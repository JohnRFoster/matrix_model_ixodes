##---------------------------------------------------------------------------------##
##                           Daily Cary tick data for JAGS                         ##
##---------------------------------------------------------------------------------##

#' Daily Cary tick data for JAGS / NIMBLE
#'
#' This function reads the raw tick drag data from the Cary institute and formats it to be
#' ready for JAGS. This data is for running models to estimate the state (abundance of
#' ticks) for each life stage at each site at the user specified interval, output data list also includes
#' the met array
#'
#' @export
#' @examples ## to estimate every week
#' @examples ## sites = c("Green Control","Henry Control","Tea Control")
#' @examples cary_ticks_JAGS(sites = sites,state.null = 7)


cary_ticks_met_JAGS <- function(sites = c("Green Control","Henry Control","Tea Control"),
                                state.interval=NULL){
  
  N_site <- length(sites)               # number of sites
  raw.dat <- read.csv("Data/Cary_tick.csv")   # read in data  
  met <- read.csv("Data/Cary_weather.csv")
  
  raw.dat$DATE <- as.Date(raw.dat$DATE) # convert to date
  
  # storage
  date.seq <- list(length(sites)) # all days within sampling window
  N_days <- vector()              # total number of days within sampling window
  df <- list()                    # number of days between sampling events
  dt.index <- list()              # number of days elaplsed from the first sampling occasion
  tick <- list(length(sites))     # data
  samp.dates <- list()            # sampling days
  interval.seq <- list()          # days to estimate latent state
  N_est <- vector()               # total number of days to estimate the latent state
  
  # separate by site
  for(j in 1:length(sites)){
    t <- subset(raw.dat, Grid == sites[j])
    t <- t[,c("DATE", "n_larvae", "n_nymphs", "n_adults")]
    tick[[j]] <- t
    samp.dates[[j]] <- t$DATE
    date.seq[[j]] <- data.frame(seq.Date(t$DATE[1], t$DATE[nrow(t)], 1))
    N_days[j] <- nrow(date.seq[[j]])
    
    # if we want to estimate the state between sampling occasions
    if(!is.null(state.interval)){
      interval <- seq(1,by=state.interval,length=floor(nrow(date.seq[[j]]) / state.interval))
      interval.seq[[j]] <- sort(as.Date(c(samp.dates[[j]],date.seq[[j]][interval,]),format = "%Y-%m-%d"))
      interval.seq[[j]] <- unique(interval.seq[[j]])
      N_est[j] <- length(interval.seq[[j]])
      df[[j]] <- as.numeric(diff(interval.seq[[j]]))
      index <- c(df[[j]],0)
      dt.index[[j]] <- cumsum(index)
      
      # if we want to estimate the state on just the sampling days
    } else {
      N_est[j] <- nrow(t)
      interval.seq[[j]] <- samp.dates[[j]]
      df[[j]] <- as.numeric(diff(t$DATE))
      index <- c(df[[j]],0)
      dt.index[[j]] <- cumsum(index)
    }
  }
  
  df.mat <- matrix(NA,length(sites),max(sapply(df,length)))
  for(j in 1:length(sites)){
    for(t in 1:length(df[[j]]))
      df.mat[j,t] <- df[[j]][[t]]
  }
  dt.index.mat <- matrix(NA,length(sites),max(sapply(dt.index,length)))
  for(j in 1:length(sites)){
    for(t in 1:length(dt.index[[j]]))
      dt.index.mat[j,t] <- dt.index[[j]][[t]]
  }
  
  # make indexing work for combining projection matrices
  df.mat[,1] <- df.mat[,1] - 1
  dt.index.mat <- cbind(rep(1, 3), dt.index.mat)
  
  # all tick array
  # dim 1 = life stage (larvae, nymph, adult)
  # dim 2 = date
  # dim 3 = site
  all.tick <- array(NA,dim = c(3,max(N_est),N_site))
  
  # match tick date to interval.seq, place in array
  for(s in 1:N_site){
    for(i in 1:nrow(tick[[s]])){
      larv <- tick[[s]][i,"n_larvae"]
      nymph <- tick[[s]][i,"n_nymphs"]
      adult <- tick[[s]][i,"n_adults"]
      col.index <- which(tick[[s]]$DATE[i] == interval.seq[[s]])
      all.tick[1,col.index,s] <- larv
      all.tick[2,col.index,s] <- nymph
      all.tick[3,col.index,s] <- adult
    }
  }
  
  
  met <- met[, c("DATE","MIN_TEMP","MAX_TEMP","MAX_RH","MIN_RH","TOT_PREC")]
  
  ## cumulative gdd calculation
  met$DATE <- lubridate::mdy(met$DATE)
  met <- met %>% 
    filter(DATE >= "1995-01-01") %>% 
    mutate(year = year(DATE))
  
  year <- met %>% 
    pull(year) %>% 
    unique()
  
  # met$year <- as.numeric(as.factor(met$year))
  
  calc_gdd <- function(base=10, max, min){
    gdd <- max(mean(max, min) - base, 0)
  }
  
  gdd.vec <- vector()
  cum.gdd <- vector()
  for(i in year){
    subset <- subset(met, year == i)
    min.missing <- which(is.na(subset$MIN_TEMP))
    max.missing <- which(is.na(subset$MAX_TEMP))
    subset$MIN_TEMP[min.missing] <- mean(subset$MIN_TEMP, na.rm = TRUE)
    subset$MAX_TEMP[max.missing] <- mean(subset$MAX_TEMP, na.rm = TRUE)
    gdd <- vector()
    for(t in 1:nrow(subset)){
      gdd[t] <- calc_gdd(10, subset$MAX_TEMP[t], subset$MIN_TEMP[t]) 
    }
    cum.gdd <- c(cum.gdd, cumsum(gdd))
    gdd.vec <- c(gdd.vec, gdd)
  }
  # names(cum.gdd) <- names(year) <- met$DATE
  
  gdd <- matrix(NA, N_site, 3767)
  met$MIN_TEMP <- scale(met$MIN_TEMP, scale = TRUE)
  met$MAX_TEMP <- scale(met$MAX_TEMP, scale = TRUE)
  met$MIN_RH <- scale(met$MIN_RH, scale = TRUE)
  met$MAX_RH <- scale(met$MAX_RH, scale = TRUE)
  met$TOT_PREC <- scale(met$TOT_PREC, scale = TRUE)
  met$DATE <- as.Date(met$DATE)
  met.seq <- list()
  for(j in 1:length(date.seq)){
    start <- as.Date(date.seq[[j]][1,1])
    start <- which(met$DATE == start)
    end <- as.Date(date.seq[[j]][nrow(date.seq[[j]]),1])
    end <- which(met$DATE == end)
    met.seq[[j]] <- start:end
    gdd[j,1:length(met.seq[[j]])] <- cum.gdd[met.seq[[j]]]
  }
  
  
  
  ## met array
  # dim 1 = days
  # dim 2 = metric(1 = max temp,
  #                2 = max rh,
  #                3 = precip,
  #                4 = vpd)
  # dim 3 = site
  met.x <- array(data = NA, dim = c(max(sapply(met.seq,length)),4,3))
  for(s in 1:N_site){
    for(t in 1:length(met.seq[[s]])){
      for(i in 2:4){
        row <- met.seq[[s]][t]
        met.x[t,i-1,s] <- met[row,i]
      }
    }
  }
  
  build.met <- function(met.var, met.seq){
    met.matrix <- matrix(NA, max(sapply(met.seq,length)), 3)
    for(s in 1:N_site){
      day.site <- length(met.seq[[s]])
      met.matrix[1:day.site,s] <- met.var[met.seq[[s]]]
    }
    return(met.matrix)
  }
  
  min.temp <- build.met(met$MIN_TEMP, met.seq)
  max.temp <- build.met(met$MAX_TEMP, met.seq)
  min.rh <- build.met(met$MIN_RH, met.seq)
  max.rh <- build.met(met$MAX_RH, met.seq)
  precip <- build.met(met$TOT_PREC, met.seq)
  
  mis.min.temp <- mis.max.temp <- mis.min.rh <- mis.max.rh <- mis.vpd
  for(s in 1:N_site){
    mis.min.temp[1:length(which(is.na(met$MIN_TEMP[met.seq[[s]]]))),s] <- which(is.na(met$MIN_TEMP[met.seq[[s]]]))
    mis.max.temp[1:length(which(is.na(met$MAX_TEMP[met.seq[[s]]]))),s] <- which(is.na(met$MAX_TEMP[met.seq[[s]]]))
    mis.min.rh[1:length(which(is.na(met$MIN_RH[met.seq[[s]]]))),s] <- which(is.na(met$MIN_RH[met.seq[[s]]]))
    mis.max.rh[1:length(which(is.na(met$MAX_RH[met.seq[[s]]]))),s] <- which(is.na(met$MAX_RH[met.seq[[s]]]))
  }
  
  disgard <- function(mis.matrix){
    dis <- which(!complete.cases(mis.matrix))  
    keep <- mis.matrix[1:(dis[1]-1),]
    return(keep)
  }
  
  mis.min.temp <- disgard(mis.min.temp) 
  mis.max.temp <- disgard(mis.max.temp) 
  mis.min.rh <- disgard(mis.min.rh) 
  mis.max.rh <- disgard(mis.max.rh) 
  
  data <- list(y = all.tick,
               dt.index = dt.index.mat,
               df = df.mat,
               N_est = N_est,
               N_days = N_days,
               temp.min = min.temp,
               temp.max = max.temp,
               rh.min = min.rh,
               rh.max = max.rh,
               precip = precip,
               mis.temp.min = mis.min.temp,
               mis.temp.max = mis.max.temp,
               mis.rh.min = mis.min.rh,
               mis.rh.max = mis.max.rh,
               gdd = gdd)
  
  return(data)
  
}
