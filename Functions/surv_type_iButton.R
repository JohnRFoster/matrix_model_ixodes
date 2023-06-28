##---------------------------------------------------------------------------------##
##                            Tick Survival Data intake                            ##
##---------------------------------------------------------------------------------##

#' Tick Survival Data Intake function
#'
#' This funciton is a wrapper around the surv_type function that also reads in
#' NOAA data and reads the raw csv file containing tick survival data
#' from the Cary institute survival experiments (soil cores) and transforms
#' the data to be fed to JAGS
#'
#' @param file Path to CSV
#' @param subset one of "Fed Larvae", "Flat Larvae", "Fed Nymph", "Flat Nymph", "Overwintering Nymph"
#' @param season default = NULL, must set to "winter or "summer for larvae
#' @export
#' @examples surv.type(file, "Fed Nymph")

library(tidyverse)
library(plantecophys)
# source("Functions/vpd.R")

surv_type_iButton <- function(subset, position){
  
  met <- read.csv("../Data/TickClimate_iButtonData.csv", stringsAsFactors = FALSE)

  wrong.fd <- which(met$X...Site == "FD ")
  if(length(wrong.fd) > 0) met$X...Site[wrong.fd] <- "FD"
  met <- met %>% 
    dplyr::filter(Position == position) %>% 
    mutate(Temp_C = (Temp-32)*5/9)
  if(position == "Above"){
    met <- met %>% 
      mutate(vpd = RHtoVPD(RH_SatDrift, Temp_C))
  }
  range(met$vpd, na.rm=TRUE)
  
  # convert from factor to date
  met$Date <- as.Date(as.character(met$Date), "%m/%d/%Y")
  
  rawdata <- read.csv("../Data/Cary_Larvae_Nymph_Surv.csv")
  dat <- rawdata
  dat$Bin <- 1:nrow(dat) # create column that is just row numbers
  
  dat$Deployment_Name <- as.character(dat$Deployment_Name)
  
  # convert from factor to date
  dat$Date_Deployed <- as.Date(as.character(dat$Date_Deployed), "%m/%d/%Y")
  dat$Date_Retrieved <- as.Date(as.character(dat$Date_Retrieved), "%m/%d/%Y")
  
  dat <- subset(dat, !is.na(dat$N_Recovered))   # subset ticks with conflicting info in N_Recovered
  dat <- subset(dat, Deployment_Name == subset) # subset data to use
  
  problem.cores <-  c("B3", "H9", "S3", "J13", "D7", "M16", "Q10", "J14")
  dat <- subset(dat, !(TC_ID %in% problem.cores))
  
  start.date <- dat$Date_Deployed   # extract start date
  end.date <- dat$Date_Retrieved    # extract end date
  
  # sequence of all days within range of subsetted data
  days.seq <- seq.Date(min(start.date),
                       max(end.date), 1)
  
  temp.below <- rh.below <- vpd.below <- matrix(NA, length(days.seq), nrow(dat))
  
  ## data specs
  for(i in 1:nrow(dat)){
    dat.site <- as.character(dat[i,"Site"])
    met.sub <- met %>% 
      dplyr::filter(Site == dat.site)
    # dplyr::filter(X...Site == dat.site)
    
    temp <- met.sub %>% 
      dplyr::group_by(Date) %>% 
      dplyr::summarise(temp = mean(Temp_C, na.rm = TRUE))
    
    rh <- met.sub %>% 
      dplyr::group_by(Date) %>%
      dplyr::summarise(rh = mean(RH_SatDrift, na.rm = TRUE))
    
    vpd <- met.sub %>% 
      dplyr::group_by(Date) %>%
      dplyr::summarise(vpd = mean(vpd, na.rm = TRUE))
    
    core.start <- which(start.date[i] == purrr::as_vector(temp[,1]))
    core.end <- which(end.date[i] == purrr::as_vector(temp[,1]))
    
    temp.core <- temp[core.start:core.end,2]
    temp.core <- purrr::as_vector(temp.core)
    rh.core <- rh[core.start:core.end,2]
    rh.core <- purrr::as_vector(rh.core)
    vpd.core <- vpd[core.start:core.end,2]
    vpd.core <- purrr::as_vector(vpd.core)
    
    temp.below[1:length(temp.core),i] <- temp.core
    rh.below[1:length(temp.core),i] <- rh.core
    vpd.below[1:length(temp.core),i] <- vpd.core
  }
  
  temp.below <- apply(temp.below, 2, scale, scale = FALSE)
  rh.below <- apply(rh.below, 2, scale, scale = FALSE)
  vpd.below <- apply(vpd.below, 2, scale, scale = FALSE)
  
  # Number of days in soil
  N_Days <- as.numeric(difftime(dat$Date_Retrieved, dat$Date_Deployed))
  
  # JAGS data for fed ticks (includes number successfully molted)
  if(subset == "Flat Larvae" ||
     subset == "Flat Nymph" ||
     subset == "Overwintering Nymph"){
    y = cbind(dat$N_Deployed, dat$N_Recovered)
  } else {
    y = cbind(dat$N_Deployed, dat$N_Survived)
  }
  
  data <- list(y = y,
               N_days = N_Days,
               # site.index = as.numeric(dat$Site),
               N = nrow(y),
               # sites = unique(as.numeric(dat$Site)),
               temp = temp.below,
               rh = rh.below,
               vpd = vpd.below)
  return(data)
}