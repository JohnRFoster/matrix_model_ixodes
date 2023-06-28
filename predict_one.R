#' This is a one step ahead prediction for tick models with dormant life stages
#' 
#' 
#'
#' @param params matrix of parameters
#' @param ticks matrix of latent state estimated from nimble
#' @param driver vector of met driver in process model
#' @param data model data fed to NIMBLE
#' @param q.index the row indices for questing ticks (rows in data$y data)

library(LaplacesDemon)

predict_one <- function(params, ticks, driver, data, mice.site = NA, mna = FALSE, mice.ua = TRUE, met.site){
  
  Nmc <- nrow(params)          # number of samples
  N_est <- data$N_est          # number of observation days
  N_days <- data$N_day         # number of total days
  dt.index <- data$dt.index    # index of when observations occur in N_days
  df <- data$df                # number of days between samples occasions
  gdd <- data$gdd[1:N_days]    # cumulative gdd vector
  seq.days <- data$seq.days    # matrix of which days to permute given days
  n.ls <- nrow(data$y)         # number of life stages
  params <- as_tibble(params)
  q.index <- which(!is.na(data$y[,1]))
  d.index <- which(is.na(data$y[,1]))
  
  use.proc <- if_else(all(params$`sig[1]` == 0), FALSE, TRUE)
  
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
      if(mice.ua){
        mice.sd <- 1/sqrt(mice$mice.prec)  
      } else {
        mice.sd <- rep(0, length(mice.mean))
      }  
    }
  }
  
  # storage
  P <- array(0, dim = c(n.ls, n.ls, N_days))
  A <- array(0, dim = c(n.ls, n.ls, N_days, Nmc))
  pred <- array(dim = c(n.ls, N_est, 2, Nmc)) #(life stage, time, latent/observation, ensemble)
  for(i in 1:n.ls) pred[i,1,1,] <- ticks %>% pull(paste0("x[", i, ", 1]"))
  
  # build daily transition matrices, vectorized over ensemble members
  for(t in seq_len(N_days)){
    
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
    
    # draw transition matrix
    if(n.ls == 4){
      l2n.quest <- ifelse(x.gdd >= rho.n & x.gdd <= 2500, 1, 0) 
      A[1,1,t,] <- phi.l * (1 - l2n)
      A[1,4,t,] <- lambda
      A[2,1,t,] <- phi.l * l2n
      A[2,2,t,] <- 1 - l2n.quest
      A[3,2,t,] <- l2n.quest
      A[3,3,t,] <- phi.n * (1 - theta.n2a)
      A[4,3,t,] <- phi.n * theta.n2a
      A[4,4,t,] <- phi.a  
    } else {
      l2n.quest <- ifelse(x.gdd >= rho.n & x.gdd <= 2500, l2n, 0) 
      A[1,1,t] <- phi.l * (1 - l2n.quest)
      A[1,3,t] <- lambda[t]
      A[2,1,t] <- phi.l * l2n.quest
      A[2,2,t] <- phi.n * (1 - theta.n2a[t])
      A[3,2,t] <- phi.n * theta.n2a[t]
      A[3,3,t] <- phi.a  
    }
    
  }
  
  for(m in 1:Nmc){
    # if(Nmc %% m == 0) message(paste(m, "/", Nmc, "simulations complete.", m/Nmc*100, "%"))
    ## aggrigate transition matricies
    for(t in 1:(N_est-1)){    # loop over the number of sampling days - 1
      
      for(day in seq.days[t,1]){
        P[,,day] <- A[,,day,m] 
      } 
      
      for(day in seq.days[t,2:df[t]]){
        P[,,day] <- P[,,day+1] %*% A[,,day,m]
      }
      
      ## predict new
      TRANS <- P[,,day] # transition matrix
      
      obs <- rep(0, n.ls) # latent state
      for(ls in 1:n.ls) obs[ls] <- ticks %>% 
        slice(m) %>% 
        pull(paste0("x[", ls, ", ", t, "]"))
      
      Ex <- TRANS %*% obs
      
      if (use.proc) { # add process error?
        OMEGA <- matrix(0, n.ls, n.ls)
        OMEGA[1, 1] <- params$`sig[1]`[m]
        OMEGA[2, 2] <- params$`sig[2]`[m]
        OMEGA[3, 3] <- params$`sig[3]`[m]
        if(n.ls == 4){
          OMEGA[4, 4] <- params$`sig[4]`[m] 
        }
        
        est <- nimble::rmnorm_chol(mean = Ex, cholesky = chol(OMEGA), prec_param = 0)
      } else {
        est <- as.vector(Ex)
      }
      
      # store latent state simulations
      pred[, t+1, 1, m] <- pmax(est, 0)
      
      # latent state simulations + observation error
      pred[q.index,t+1,2,m] <- rpois(length(q.index), pmax(est[q.index], 0.01)) 
      
      if(n.ls == 3){
        pred[d.index,t+1,2,m] <- max(0, est[d.index]) # no observation model on dormant states
      }
      
    }
  }
  
  return(pred)
}


