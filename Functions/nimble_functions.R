#' build data, constants, and met to feed into nimbleModel
#' 
#' @param site.run the field site being modelled
#' 

nimble_data <- function(site.run){
  library(lubridate)
  library(tidyverse
          )
  raw.dat <- read.csv("Data/Cary_ticks.csv")   # read in data  
  
  df <- raw.dat %>% 
    as_tibble() %>% 
    filter(Grid == site.run) %>% 
    select(DATE, n_larvae, n_nymphs, n_adults) %>% 
    mutate(DATE = ymd(DATE))
  
  y <- tibble(larvae = df$n_larvae,
              d.nymph = NA,
              q.nymph = df$n_nymphs,
              adult = df$n_adults) %>% 
    t()
  
  all.days <- seq.Date(first(df$DATE), last(df$DATE), by = 1)
  N_est <- ncol(y)
  N_days <- length(all.days)
  df <- as.numeric(diff(df$DATE))
  dt.index <- c(1, cumsum(c(df, 0)))
  
  met <- read.csv("Data/Cary_weather.csv")
  met.df <- met %>% 
    as_tibble() %>% 
    mutate(DATE = mdy(DATE)) %>% 
    filter(DATE %in% all.days) %>% 
    select(DATE, MAX_TEMP, MAX_RH, MIN_RH, TOT_PREC) 
  
  max.temp <- pull(met.df, MAX_TEMP)
  xout <- which(is.na(max.temp))
  fill.mis.temp <- approx(max.temp, xout = xout)$y
  max.temp[xout] <- fill.mis.temp
  
  met.df <- met.df %>% 
    mutate(year = year(DATE),
           max.temp.filled = max.temp) %>% 
    group_by(year) %>% 
    mutate(gd = if_else(max.temp.filled > 10, max.temp.filled-10, 0),
           cumGD = cumsum(gd)) %>% 
    select(-max.temp.filled) %>% 
    mutate(MAX_TEMP = scale(MAX_TEMP)[,1],
           MAX_RH = scale(MAX_RH)[,1],
           MIN_RH = scale(MIN_RH)[,1],
           TOT_PREC = scale(TOT_PREC)[,1]) %>% 
    ungroup()
  
  seq.days <- matrix(NA, N_est - 1, max(df, na.rm = TRUE))
  for (i in 1:(N_est - 1)) {
    xx <- (dt.index[i + 1] - 1):dt.index[i]
    seq.days[i, 1:length(xx)] <- xx
  }
  
  
  data <- list()
  data$y <- y
  data$gdd <- pull(met.df, cumGD)
  
  constants <- list()
  constants$N_est <- N_est
  constants$N_days <- N_days
  df[1] <- df[1] - 1
  constants$df <- df
  constants$dt.index <- dt.index
  constants$seq.days <- seq.days
  
  return(list(data = data,
              constants = constants,
              met = met.df))
  
}


#' if_else statement for using in a nimble model
#' 
#' @param condition statement to evaluate T/F
#' @param valueIf outcome if TRUE
#' @param valueElse outcome if FALSE
#' 
if_else_nimble <- nimbleFunction(
  run = function(condition = integer(0), valueIf=double(0), valueElse=double(0)) {
    returnType(double(0))
    if (condition==TRUE) {
      return(valueIf)
    } else {
      return(valueElse)
    }
  }
)
# if_else_nimble <- compileNimble(if_else_nimble)


#' function for running nimble chains in parallel
#' this function checks convergence and will keep sampling if not converged
#' 
#' @param cl cluster built in main script
#' @param model nimble model code
#' @param constants nimble constants list
#' @param data nimble data list
#' @param inits initial values function
#' @param monitor vector of variables to monitor
#' @param file.name full file path to where you want intermediate output to be saved. If NA nothing saves.
#' @param use.dzip using zero-inflated Poisson?
#' @param n.ens the number of total iterations to save
#' @param check.interval before convergence, how often do you want to print/save output? 
#' @param thin thinning interval
#' @param n.iter the number of iterations for each intermediate sampling step
#' @param n.burnin the number of burnin iterations
#' @param psrf.max the maximum value allowed for convergence (R_hat threshold)
#' @param min.effective.size the minimum number of effective samples allowed
#' @param max.iter the maximum number of iterations before sampling stops
#' @param save.states are we going to save latent state samples? default TRUE

run_nimble_parallel <- function(cl, model, constants, data, inits, monitor, file.name = NA,
                                use.dzip = FALSE, n.ens = 5000, check.interval = 10, thin = 1, mna = NA,
                                n.iter = 5e6, n.burnin = 0, psrf.max = 1.1, min.effective.size = 3000,
                                max.iter = 3e6, calculate = TRUE, use.conjugacy = TRUE, rjmcmc = FALSE){
  library(parallel)
  library(nimble)
  library(coda)
  
  n.cores <- length(cl) # number of cores used
  
  
  export.vec <- c(
    "model",
    "constants",
    "data",
    "n.iter",
    "n.burnin",
    "monitor",
    "thin",
    "n.cores",
    "if_else_nimble",
    "calculate",
    "use.conjugacy",
    "rjmcmc",
    "mna"
  )
  
  # export everything we need to cluster nodes
  if(use.dzip){
    source("Functions/ZIP.R")
    assign("dZIP", dZIP, envir = .GlobalEnv)
    assign("rZIP", rZIP, envir = .GlobalEnv)
    export.vec <- c(export.vec, "dZIP", "rZIP")
  }
  
  clusterExport(cl,
                export.vec,
                envir = environment()) 
  
  
  # export inits to clusters
  for(j in seq_along(cl)){
    set.seed(j)
    init <- inits()
    # print(init)
    clusterExport(cl[j], "init", envir = environment())
  }
  
  message("Running mcmc...")
  start.time <- Sys.time()
  out <- clusterEvalQ(cl, { # sample on each cluster
    library(nimble)
    library(nimbleEcology)
    library(coda)
    model.rw <- nimbleModel(model,
                            constants = constants,
                            data = data,
                            inits = init,
                            calculate = calculate)
    cModel.rw <- compileNimble(model.rw)
    mcmcConf <- configureMCMC(cModel.rw, 
                              monitors = monitor,
                              thin = thin,
                              useConjugacy = use.conjugacy)
    
    mcmcConf$removeSampler(target = c(
      "phi.a.mu",
      "phi.n.mu",
      "phi.l.mu",
      "theta.na",
      "theta.ln",
      "sig"
    ))
    mcmcConf$addSampler(target = "phi.a.mu", type = "ess")
    mcmcConf$addSampler(target = "phi.n.mu", type = "ess")
    mcmcConf$addSampler(target = "phi.l.mu", type = "ess")
    mcmcConf$addSampler(target = "theta.ln", type = "ess")
    mcmcConf$addSampler(target = "theta.na", type = "ess")
    mcmcConf$addSampler(target = "sig", type = "AF_slice")
    
    if(rjmcmc){
      configureRJ(mcmcConf,
                  targetNodes = "beta",
                  indicatorNodes = "z")
    }
    mcmcBuild <- buildMCMC(mcmcConf)
    compMCMC <- compileNimble(mcmcBuild)
    compMCMC$run(niter = n.iter, nburnin = n.burnin)
    return(as.mcmc(as.matrix(compMCMC$mvSamples)))
  })
  delta.time <- Sys.time() - start.time
  message(paste(n.iter,"iterations in"))
  print(delta.time)
  out.mcmc <- as.mcmc.list(out)
  check.mcmc <- out.mcmc
  
  states.monitor <- if_else(any(grepl("x[", colnames(out.mcmc[[1]]), fixed = TRUE)),
                            TRUE, FALSE)
  z.monitor <- if_else(any(grepl("z[", colnames(out.mcmc[[1]]), fixed = TRUE)),
                       TRUE, FALSE)
  
  if(states.monitor){
    check.mcmc <- check.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
  }
  
  # check for convergence and effective sample size
  conv_eff <- function(xx, z){
    g.diag <- gelman.diag(xx, multivariate = FALSE)$psrf
    effect.sizes <- effectiveSize(xx)
    if(z){
      g.diag <- g.diag[-grep("z[", rownames(g.diag), fixed = TRUE),]
      effect.sizes <- effect.sizes[-grep("z[", names(effect.sizes), fixed = TRUE)]
    }
    return(list(g.diag = g.diag, effect.sizes = effect.sizes))
  }
  
  mcmc.check <- conv_eff(check.mcmc, z.monitor)
  g.diag <- mcmc.check$g.diag
  effect.sizes <- mcmc.check$effect.sizes
  
  convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
  message(paste("Convergence:", convergence)) 
  if(is.na(convergence)){
    print(g.diag)
    convergence <- FALSE
  } 
  
  enough.samples <- min(effect.sizes) >= min.effective.size
  message(paste("Enough Effective Samples:", enough.samples)) 
  
  counter <- 1
  total.iter <- counter * n.iter
  resetMV <- FALSE
  while(!convergence | !enough.samples){
    counter <- counter + 1
    total.iter <- counter * n.iter
    
    message("Running mcmc...")
    start.time <- Sys.time()
    clusterExport(cl, "resetMV", envir = environment())
    out2 <- clusterEvalQ(cl, {
      compMCMC$run(n.iter, reset = FALSE, resetMV = resetMV)
      return(as.mcmc(as.matrix(compMCMC$mvSamples)))
    })
    out.mcmc <- as.mcmc.list(out2)
    
    delta.time <- Sys.time() - start.time
    message(paste(n.iter,"iterations in"))
    print(delta.time)
    
    if(states.monitor){
      check.mcmc <- out.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
    } else {
      check.mcmc <- out.mcmc
    }
    
    mcmc.check <- conv_eff(check.mcmc, z.monitor)
    g.diag <- mcmc.check$g.diag
    convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
    if(is.na(convergence)) convergence <- FALSE
    message(paste("Convergence:", convergence))  
    
    effect.sizes <- mcmc.check$effect.sizes
    enough.samples <- min(effect.sizes) >= min.effective.size
    message(paste("Enough Effective Samples:", enough.samples))
    
    message(paste("Total iterations:", total.iter))
    if(total.iter > max.iter) break
    
    if(counter %% check.interval == 0) {
      message("Convergence check:")
      print(g.diag)
      message("Effective sample sizes:")
      print(effect.sizes)
      if(!is.na(file.name)){
        
        if(convergence){
          GBR <- gelman.plot(check.mcmc, ask = FALSE)
          burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
          message(paste("Burnin after", burnin, "iterations"))
          out.mcmc <- window(as.mcmc.list(out.mcmc), start = burnin)
        }
        
        ## split output
        mfit <- as.matrix(out.mcmc, chains = TRUE)
        pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
        chain.col <- which(colnames(mfit) == "CHAIN")
        predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
        params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
        
        save.ls <- list(params = params,
                        predict = predict,
                        convergence = convergence,
                        effect.sizes = effect.sizes,
                        model.code = model.code)
        
        # if(max(g.diag[,"Point est."]) < psrf.max){
        #   save.ls <- purrr::prepend(
        #     save.ls,
        #     list(predict = predict))
        # }
        
        message("Writing output...")
        save(save.ls, file = file.name)
        message("DONE")
      }
      
      resetMV <- if_else(max(g.diag[, "Upper C.I."]) > 20,
                         TRUE,
                         FALSE)
    }
    
  }
  
  if(convergence & enough.samples){ #
    GBR <- gelman.plot(check.mcmc, ask = FALSE)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
    message(paste("Burnin after", burnin, "iterations"))
    
    # sometimes burnin happens at the end of the chain, need to keep sampling if that happens
    if(total.iter - burnin < 1000){
      message("Additional iterations because burnin is too close to the end of the chain")
      out3 <- clusterEvalQ(cl, {
        compMCMC$run(n.iter, reset = FALSE, resetMV = TRUE)
        return(as.mcmc(as.matrix(compMCMC$mvSamples)))
      })
      out.mcmc <- as.mcmc.list(out.3)
    } else{
      out.mcmc <- window(as.mcmc.list(out.mcmc), start = burnin)
    }  
    
    # return a thinned matrix (raw mcmc objects can get big)
    samples <- as.matrix(out.mcmc)
    if(n.ens < nrow(samples)){
      thin.seq <- round(seq(1, nrow(samples), length.out = n.ens)) # sequence of samples to keep
      samples <- samples[thin.seq,]  
    }
    
  } else {
    ## split output
    mfit <- as.matrix(out.mcmc, chains = TRUE)
    pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
    chain.col <- which(colnames(mfit) == "CHAIN")
    predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
    params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
    samples <- list()
    samples$predict <- predict
    samples$params <- params
  }
  
  return(list(samples = samples,
              convergence = convergence,
              enough.samples = enough.samples)) 
  
}

