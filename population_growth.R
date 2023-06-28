


# from table 14.2 in caswell

#' calculate population growth rate
#'
#' @param matrices the daily transition matix
#' 
mcmc_lam <- function(matrices){
  if (is.list(matrices)) {
    matrices <- matrix(unlist(matrices), ncol = length(matrices))
  }
  s <- sqrt(dim(matrices)[1])
  n.days <- ncol(matrices)
  n0 <- rep(1, s) / s
  r <- numeric(n.days)
  for(t in 1:n.days){
    A <- matrix(matrices[, t], nrow = s)
    n0 <- A %*% n0
    N <- sum(n0)
    r[t] <- log(N)
    n0 <- n0/N
  }
  mean(r)  
}

ticks <- read_csv("Data/Cary_ticks.csv")
dates <- ticks %>% 
  filter(Grid == site.run) %>% 
  mutate(year = year(DATE)) %>% 
  group_by(year) %>% 
  mutate(seq = n():1) %>% 
  select(DATE, seq)

#' calculate population growth rate over all periods and mcmc samples
#'
#' @param matrices the daily transition matix
#' 
simulate_pop_growth <- function(A.ls, data){
  A <- A.ls$A
  Nmc <- dim(A)[4]
  N_est <- data$N_est
  N_days <- data$N_day
  dt.index <- data$dt.index
  df <- data$df
  gdd <- data$gdd[1:N_days]
  n.ls <- 4
  seq.days <- data$seq.days
  samp.date <- dates$DATE[-1]
  pop.metrics <- tibble()
  for(m in 1:Nmc){
    for(t in 1:(N_est-1)){   # loop over the number of sampling days - 1
      x <- 1
      S <- matrix(NA, 16, df[t])
      index <- df[t]:1
      for(day in seq.days[t,1:(df[t])]){
        S[,index[x]] <- as.vector(A[,,day,m])
        x <- x + 1
      }
      stoch.sim <- mcmc_lam(S)
      
      lambdas <- tibble(
        time = samp.date[t],
        mcmc = m,
        lambda.sim = exp(stoch.sim),
        stage = "population"
      )
      pop.metrics <- bind_rows(pop.metrics, lambdas)
    }
    if (m == 1 || m %% 100 == 0) {
      message("  mcmc ", m)
    }
  }
  return(pop.metrics)
}

