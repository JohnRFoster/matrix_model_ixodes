#' get mouse abundance estimates from jags output
#'
#' @param site.run filed site mice were estimated

mice_estimated_jags <- function(site.run){
  
  ### mouse data
  mice <- read.csv("Data/Cary_mouse.csv")
  mice <- subset(mice, Grid == site.run)
  mice$Full.Date.1 <- as.Date(mice$Full.Date.1)
  mice$Full.Date.2 <- as.Date(mice$Full.Date.2)
  
  # unique trap 
  day.1 <- as.character(unique(mice$Full.Date.1)) # 1st capture event in mouse series
  day.2 <- as.character(unique(mice$Full.Date.2)) # 2nd capture event in mouse series
  
  # get dates
  day.mice <- c(rbind(day.1, day.2)) 
  mice.obs <- as.Date(day.mice) # unique sampling days: mice
  
  # every day in mouse sequence
  mice.seq <- seq.Date(mice.obs[1], mice.obs[length(mice.obs)], by = 1)
  
  # load mice model output
  if(site.run == "Green Control"){
    load("FinalOut/GreenControlMR/GreenControlMR_1_5.RData")
  } else if(site.run == "Henry Control"){
    load("FinalOut/HenryControlMR/HenryControlMR_1_6.RData")
  } else {
    load("FinalOut/TeaControlMR/TeaControlMR_1_6.RData")
  }
  
  predict <- as.matrix(out$predict)
  
  # thin to 5000 iterations
  predict <- predict[seq(1, nrow(predict), by = round(nrow(predict)/5000)),] 
  # predict <- apply(predict, 2, sort)
  
  mice.est.ci <- list()
  mice.est.ci[[1]] <- predict[,1]
  for(t in 2:ncol(predict)){
    time.diff <- as.numeric(difftime(mice.obs[t], mice.obs[t-1]))
    if(time.diff > 1){
      x <- c(1, time.diff)
      interpol <- matrix(NA, nrow(predict), time.diff+1)
      for(i in 1:nrow(predict)){
        y <- c(predict[i,t-1], predict[i,t])
        interpol[i,] <- approx(x, y, n = time.diff+1, method = "linear")$y 
      }
      mice.est.ci[[t]] <- interpol 
    } 
  }
  
  # create data frame and calculate mean and standard dev
  mice.sort.all.days <- do.call(cbind, mice.est.ci)
  
  mice.int <- round(mice.sort.all.days, 0)
  
  mice.all.days.mean <- apply(mice.sort.all.days, 2, mean)
  mice.all.days.var <- apply(mice.sort.all.days, 2, var)
  
  # tick data
  dat <- read.csv("Data/Cary_ticks.csv") # tick data
  tick <- dat[,c("Grid", "DATE", "n_larvae", "n_nymphs", "n_adults")]
  tick <- subset(tick, Grid == site.run)
  tick$DATE <- as.Date(tick$DATE)
  
  # match tick dates to mice dates
  start.tick <- which(tick$DATE[1] == mice.seq)
  end.tick <- which(tick$DATE[length(tick$DATE)] == mice.seq)
  
  # index mice estimates - convert sd to prec
  mice.for.jags.mean <- mice.all.days.mean[start.tick:end.tick]
  mice.for.jags.var <- mice.all.days.var[start.tick:end.tick]
  
  mice.scale.mean <- scale(mice.for.jags.mean)
  mice.scale.var <- (1/sd(mice.for.jags.mean))^2 * mice.for.jags.var
  
  return(list(mice.mean = mice.scale.mean,
              mice.prec = 1 / mice.scale.var))
}





