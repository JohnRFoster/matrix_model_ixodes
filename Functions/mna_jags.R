#' get mouse minimum number alive (mna)
#'
#' @param site.run filed site mice were estimated


mna_jags <- function(site.run, return.mean = FALSE){
  source("Functions/ch_cary.R")
  source("Functions/known_states.R")
  if(!grepl("Control", site.run)) site.run <- paste(site.run, "Control")
  ch <- suppressWarnings(ch_cary(site.run))
  ks <- known_states(ch)
  mna <- apply(ks, 2, sum)
  
  mice.obs <- ymd(colnames(ch)) # unique sampling days: mice
  
  # every day in mouse sequence
  mice.seq <- seq.Date(mice.obs[1], mice.obs[length(mice.obs)], by = 1)
  
  mna.all.days <- rep(NA, length(mice.seq))
  mna.count <- 1
  for(i in seq_along(mice.seq)){
    if(mice.seq[i] %in% mice.obs){
      mna.all.days[i] <- mna[mna.count]
      mna.count <- mna.count + 1
    } else {
      mna.all.days[i] <- mna[mna.count]
    }
  }
  
  # tick data
  dat <- read.csv("Data/Cary_ticks.csv") # tick data
  tick <- dat[,c("Grid", "DATE", "n_larvae", "n_nymphs", "n_adults")]
  tick <- subset(tick, Grid == site.run)
  tick$DATE <- as.Date(tick$DATE)
  
  # match tick dates to mice dates
  start.tick <- which(tick$DATE[1] == mice.seq)
  end.tick <- which(tick$DATE[length(tick$DATE)] == mice.seq)
  
  # index mice estimates - convert sd to prec
  mna.for.jags <- mna.all.days[start.tick:end.tick]
  
  # center and scale
  mna.scaled <- scale(mna.for.jags)
  
  if(return.mean){
    return(
      list(
        mna = mna.scaled,
        mean = mean(mna.for.jags),
        sd = sd(mna.for.jags)
      )
    )
  } else{
    return(mna.scaled)  
  }
  
}

