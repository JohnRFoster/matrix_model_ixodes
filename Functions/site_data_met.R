##---------------------------------------------------------------------------------##
##                           Subset Daily Cary tick data for JAGS                  ##
##---------------------------------------------------------------------------------##

#' Daily Cary tick data for JAGS / NIMBLE
#'
#' This function reads takes the output from cary_tick_met_JAGS and subsets it to a specific site
#'



library(dplyr)
library(lubridate)

site_data_met <- function(site, met.variable = TRUE, data, dir = "", 
                          obs.model = FALSE, time.effect = NULL){
  # index for site
  dat.jags <- list()
  if(site == "Green Control"){
    s <- 1
    dat.jags$y <- data$y[,-73,1]
  } else if(site == "Henry Control"){
    s <- 2
    dat.jags$y <- data$y[,,2]
  } else {
    s <- 3
    dat.jags$y <- data$y[,-73,3]
  }
  
  dat.jags$N_est <- data$N_est[s]
  dat.jags$N_days <- data$N_days[s]
  dat.jags$dt.index <- data$dt.index[s,]
  dat.jags$df <- data$df[s,]
  dat.jags$gdd <- data$gdd[s,1:dat.jags$N_days]
  # dat.jags$year <- data$year[,s]
  
  if(met.variable){
    dat.jags$temp.max <- data$temp.max[1:dat.jags$N_days,s]
    dat.jags$temp.max.mis <- data$mis.temp.max[,s]
    dat.jags$temp.max.range <- range(dat.jags$temp.max, na.rm = TRUE)
    
    dat.jags$temp.min <- data$temp.min[1:dat.jags$N_days,s]
    dat.jags$temp.min.mis <- data$mis.temp.min[,s]
    dat.jags$temp.min.range <- range(dat.jags$temp.min, na.rm = TRUE)
    
    dat.jags$rh.max <- data$rh.max[1:dat.jags$N_days,s]
    dat.jags$rh.max.mis <- data$mis.rh.max[,s]
    dat.jags$rh.max.range <- range(dat.jags$rh.max, na.rm = TRUE)
    
    dat.jags$rh.min <- data$rh.min[1:dat.jags$N_days,s]
    dat.jags$rh.min.mis <- data$mis.rh.min[,s]
    dat.jags$rh.min.range <- range(dat.jags$rh.min, na.rm = TRUE)
    
    dat.jags$precip <- data$precip[1:dat.jags$N_days,s]
  }
  
  if(!is.null(time.effect)){
    raw.dat <- read.csv("Data/tick_cleaned")
    raw.dat$DATE <- as.Date(raw.dat$DATE) # convert to date
    
    dates <- raw.dat %>% 
      filter(Grid == site) %>% 
      select(DATE)
    if(time.effect == "month"){
      dat.jags$month.index <- lubridate::month(dates[,1])
      dat.jags$n.months <- unique(lubridate::month(dates[,1]))
    }  
    if(time.effect == "year"){
      dat.jags$year.index <- lubridate::year(dates[,1]) - 1993
      dat.jags$n.years <- unique(lubridate::year(dates[,1])) - 1993
    }
  }
  
  if(obs.model){
    obs.index <- c(1, dat.jags$dt.index)
    met.obs <- data$temp.min[obs.index[1:dat.jags$N_est],s]
    met.obs.miss <- which(is.na(met.obs))
    met.obs.range <- range(met.obs, na.rm = TRUE)
    
    dat.jags$met.obs <- met.obs
    dat.jags$met.obs.miss <- met.obs.miss
    dat.jags$met.obs.range <- met.obs.range
  }
  
  seq.days <- matrix(NA, dat.jags$N_est - 1, max(dat.jags$df, na.rm = TRUE))
  for (i in 1:(dat.jags$N_est - 1)) {
    xx <- (dat.jags$dt.index[i + 1] - 1):dat.jags$dt.index[i]
    seq.days[i, 1:length(xx)] <- xx
  }
  dat.jags$seq.days <- seq.days
  
  # dat.jags$R <- diag(1, 3, 3)
  
  return(dat.jags)
}