#' create ensembles of daily transition matrices for population growth rate analysis
#'
#' @param params matrix of parameter samples
#' @param driver matrix of weather drivers
#' @param data data list fed to nimble
#' @param mice.site the site mice were estimated
#' @param mna was the mice data assimilated mna or estimated from nimble?
#' @param met.site the site the model was fit

A_dormant <- function(params, driver, data, mice.site = NA, mna = FALSE, met.site){
  
  Nmc <- nrow(params)          # number of samples
  N_est <- data$N_est          # number of observation days
  N_days <- data$N_day         # number of total days
  dt.index <- data$dt.index    # index of when observations occur in N_days
  df <- data$df                # number of days between samples occasions
  gdd <- data$gdd[1:N_days]    # cumulative gdd vector
  seq.days <- data$seq.days    # matrix of which days to permute given days
  n.ls <- 4                    # number of life stages
  params <- as_tibble(params)
  
  if(!is.na(mice.site)){
    # get mice data
    if(mna){
      source("Functions/mna_jags.R")
      mice <- mna_jags(mice.site)
      mice.mean <- as.vector(mice) 
      mice.sd <- rep(0, length(mice.mean))
    } else {
      source("Functions/mice_estimated_jags.R")
      mice <- mice_estimated_jags(paste(mice.site, "Control"))
      mice.mean <- mice$mice.mean %>% as.vector()
      mice.sd <- 1/sqrt(mice$mice.prec)  
    }
  }
  
  # storage
  A <- array(0, dim = c(n.ls, n.ls, N_days, Nmc))
  phi.1 <- matrix(0, N_days, Nmc)
  phi.2 <- matrix(0, N_days, Nmc)
  phi.3 <- matrix(0, N_days, Nmc)
  theta.1 <- matrix(0, N_days, Nmc)
  theta.2 <- matrix(0, N_days, Nmc)
  theta.3 <- matrix(0, N_days, Nmc)
  lambda.m <- matrix(0, N_days, Nmc)
  
  # build daily transition matrices, vectorized over ensemble members
  for(t in seq_len(N_days)){
    
    # if(!is.na(driver)) x <- pull(driver, t)
    
    if(is.na(met.site)){
      phi.l <- invlogit(params$phi.l.mu)
      phi.n <- invlogit(params$phi.n.mu)
      phi.a <- invlogit(params$phi.a.mu)
    } else {
      xx <- slice(driver, t)
      X <- tibble(MAX_TEMP = rep(pull(xx, "MAX_TEMP"), Nmc),
                  MAX_RH = rep(pull(xx, "MAX_RH"), Nmc),
                  MIN_RH = rep(pull(xx, "MIN_RH"), Nmc),
                  TOT_PREC = rep(pull(xx, "TOT_PREC"), Nmc))
      
      if(any(is.na(xx))){
        fill.na <- which(is.na(xx))
        for(f in seq_along(fill.na)){
          X[,fill.na[f]] <- rnorm(Nmc)  
        }
      } 
      
      
      phi.l <- invlogit(params$phi.l.mu +
                          params$`beta[1]` * X$`MAX_TEMP` +
                          params$`beta[2]` * X$`MAX_RH` +
                          params$`beta[3]` * X$`MIN_RH` +
                          params$`beta[4]` * X$`TOT_PREC`)
      phi.n <- invlogit(params$phi.n.mu +
                          params$`beta[5]` * X$`MAX_TEMP` +
                          params$`beta[6]` * X$`MAX_RH` +
                          params$`beta[7]` * X$`MIN_RH` +
                          params$`beta[8]` * X$`TOT_PREC`)
      phi.a <- invlogit(params$phi.a.mu +
                          params$`beta[9]` * X$`MAX_TEMP` +
                          params$`beta[10]` * X$`MAX_RH` +
                          params$`beta[11]` * X$`MIN_RH` +
                          params$`beta[12]` * X$`TOT_PREC`)
    
    }

    rho.l <- 1400
    rho.n <- 400
    rho.a <- 2500
    
    if(!is.na(mice.site)){
      mice <- rnorm(Nmc, mice.mean[t], mice.sd[t]) 
      if("beta.m[1]" %in% colnames(params)){
        l2n <- invlogit(params$theta.ln + params$`beta.m[1]`*mice)
        n2a <- invlogit(params$theta.na + params$`beta.m[2]`*mice)  
      } else {
        l2n <- invlogit(params$theta.ln + params$`beta[13]`*mice)
        n2a <- invlogit(params$theta.na + params$`beta[14]`*mice) 
      }
    } else {
      l2n <- invlogit(params$theta.ln)
      n2a <- invlogit(params$theta.na)
    }
    
    x.gdd <- rep(data$gdd[t], Nmc)
    lambda <- ifelse(x.gdd >= rho.l & x.gdd <= 2500, params$repro.mu, 0)
    theta.n2a <- ifelse(x.gdd <= 1000 | x.gdd >= rho.a, n2a, 0)
    l2n.quest <- ifelse(x.gdd >= rho.n & x.gdd <= 2500, 1, 0)  
    
    # draw transition matrix
    A[1,1,t,] <- phi.l * (1 - l2n)
    A[1,4,t,] <- lambda
    A[2,1,t,] <- phi.l * l2n
    A[2,2,t,] <- 1 - l2n.quest
    A[3,2,t,] <- l2n.quest
    A[3,3,t,] <- phi.n * (1 - theta.n2a)
    A[4,3,t,] <- phi.n * theta.n2a
    A[4,4,t,] <- phi.a
    
    # save lower level parameters
    phi.1[t,] <- phi.l
    phi.2[t,] <- phi.n
    phi.3[t,] <- phi.a
    theta.1[t,] <- l2n
    theta.2[t,] <- l2n.quest
    theta.3[t,] <- theta.n2a
    lambda.m[t,] <- lambda
    
  }
  lower <- list(phi.1 = phi.1,
                phi.2 = phi.2,
                phi.3 = phi.3,
                theta.1 = theta.1,
                theta.2 = theta.2,
                theta.3 = theta.3,
                lambda.m = lambda.m)
  return(list(A = A,
              lower = lower))
}
