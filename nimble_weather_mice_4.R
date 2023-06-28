library(nimble)
source("Functions/nimble_functions.R")

model.code <- nimbleCode({
  
  ### priors
  phi.l.mu ~ dnorm(larva.mean, tau = larva.prec) # larvae survival
  phi.n.mu ~ dnorm(nymph.mean, tau = nymph.prec) # nymph survival
  phi.a.mu ~ dnorm(5, tau = 1)                   # adult survival
  theta.ln ~ dnorm(ln.mean, tau = 1)             # larvae -> dormant nymph daily transition 
  theta.na ~ dnorm(na.mean, tau = 1)                 # larvae -> questing nymph daily transition 
  repro.mu ~ T(dnorm(30, tau = 0.01), 0, Inf)     # reproduction
  lambda.b ~ dinvgamma(1, 1)
  
  for(j in 1:11){
    beta[j] ~ dnorm(0, tau = lambda.b)
    # beta[j] ~ dnorm(0, tau = 1/1.5^2) 
  }
  
  for(i in 1:ns){
    sig[i] ~ dinvgamma(0.01, 0.01)  
  }
  
  ### precision priors
  OMEGA[1,1] <- sig[1]
  OMEGA[1,2] <- 0
  OMEGA[1,3] <- 0
  OMEGA[1,4] <- 0
  OMEGA[2,1] <- 0
  OMEGA[2,2] <- sig[2]
  OMEGA[2,3] <- 0
  OMEGA[2,4] <- 0
  OMEGA[3,1] <- 0
  OMEGA[3,2] <- 0
  OMEGA[3,3] <- sig[3]
  OMEGA[3,4] <- 0
  OMEGA[4,1] <- 0
  OMEGA[4,2] <- 0
  OMEGA[4,3] <- 0
  OMEGA[4,4] <- sig[4]
  
  # Cholesky decomposition
  Ochol[1:ns,1:ns] <- chol(OMEGA[1:ns,1:ns])
  
  ### first latent process 
  ### first latent process 
  x[1, 1] ~ dnorm(x.init[1], tau = 0.01)
  x[2, 1] ~ dnorm(x.init[2], tau = 0.01)
  x[3, 1] ~ dnorm(x.init[3], tau = 0.01)
  x[4, 1] ~ dnorm(x.init[4], tau = 0.01)
  
  # missing X prior
  for(j in 1:n.X1){
    X[X1.mis[j],1] ~ dnorm(0, tau = 1)
  }
  for(j in 1:n.X2){
    X[X2.mis[j],2] ~ dnorm(0, tau = 1)
  }
  for(j in 1:n.X3){
    X[X3.mis[j],3] ~ dnorm(0, tau = 1)
  }
  
  ### define parameters
  for(t in 1:N_days){   # loop over every day in time series
    
    mice[t] ~ dnorm(mice.mean[t], tau = mice.prec[t]) 
    logit(l2n[t]) <- theta.ln + beta[10]*mice[t]
    logit(n2a[t]) <- theta.na + beta[11]*mice[t]
    
    theta.n2a[t] <- if_else_nimble((gdd[t] <= 1000) | (gdd[t] >= 2500),n2a[t],0)
    lambda[t] <- if_else_nimble((gdd[t] >= 1400) & (gdd[t] <= 2500),repro.mu,0)
    l2n.quest[t] <- if_else_nimble((gdd[t] >= 400) & (gdd[t] <= 2500),1,0)
    
    # daily survival model
    logit(phi.l[t]) <- phi.l.mu + 
      beta[1]*X[t,1] + # MAX_TEMP
      beta[1]*X[t,2] + # MAX_RH
      beta[2]*X[t,3] + # MIN_RH
      beta[3]*X[t,4]   # TOT_PREC
    
    logit(phi.n[t]) <- phi.n.mu + 
      beta[4]*X[t,1] + # MAX_TEMP
      beta[5]*X[t,2] + # MAX_RH
      beta[5]*X[t,3] + # MIN_RH
      beta[6]*X[t,4]   # TOT_PREC
    
    logit(phi.a[t]) <- phi.a.mu + 
      beta[7]*X[t,1] + # MAX_TEMP
      beta[8]*X[t,3] + # MIN_RH
      beta[9]*X[t,4]   # TOT_PREC
    
    A[1,1,t] <- phi.l[t] * (1 - l2n[t])
    A[1,2,t] <- 0
    A[1,3,t] <- 0
    A[1,4,t] <- lambda[t]
    A[2,1,t] <- phi.l[t] * l2n[t]
    A[2,2,t] <- 1 - l2n.quest[t]
    A[2,3,t] <- 0
    A[2,4,t] <- 0
    A[3,1,t] <- 0
    A[3,2,t] <- l2n.quest[t]
    A[3,3,t] <- phi.n[t] * (1 - theta.n2a[t])
    A[3,4,t] <- 0
    A[4,1,t] <- 0 
    A[4,2,t] <- 0
    A[4,3,t] <- phi.n[t] * theta.n2a[t]
    A[4,4,t] <- phi.a[t]
  }
  
  ### Process Model
  for(t in 1:(N_est-1)){
    
    P[1:ns,1:ns,seq.days[t,1]] <- A[1:ns,1:ns,seq.days[t,1]]
    
    for(day in 2:df[t]){
      P[1:ns,1:ns,seq.days[t,day]] <- P[1:ns,1:ns,seq.days[t,day]+1] %*% A[1:ns,1:ns,seq.days[t,day]]
    }
    
    # expected number questing
    Ex[1:ns,t] <- P[1:ns,1:ns,dt.index[t]] %*% x[1:ns,t] 
    
    # process error
    px[1:ns,t] ~ dmnorm(mean = Ex[1:ns,t], cholesky = Ochol[1:ns,1:ns], prec_param = 0)
    
    x[1,t+1] <- px[1,t]
    x[2,t+1] <- max(px[2,t], 0)
    x[3,t+1] <- px[3,t]
    x[4,t+1] <- px[4,t]
    
  }
  
  ### Data Model ###
  for(t in 1:N_est){
    
    ## fit the blended model to observed data 
    y[1,t] ~ dpois(x[1,t])
    y[3,t] ~ dpois(x[3,t])
    y[4,t] ~ dpois(x[4,t])
    
  } # t
  
})

monitor <- c("x",
             "phi.l.mu",
             "phi.n.mu",
             "phi.a.mu",
             "theta.ln",
             "theta.na",
             "repro.mu",
             "sig",
             "lambda.b"
             "beta")


