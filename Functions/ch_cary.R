##---------------------------------------------------------------------------------##
##   Mouse data intake to make CAPTURE HISTORY MATRIX from Cary Institute data     ##
##---------------------------------------------------------------------------------##

#' Capture History Matrix Function
#'
#' This function reads the raw csv file containing white footef mice capture data
#' @param path File path to csv, Default = "" and assumes csv is in working directory
#' @param grid What grid do you want to create matrix for? One of Green Control, Henry Control, Tea Control, Green Experimental, Henry Experimental, Tea Experimental
#' @examples ch.cary("Green Control")
#' @export

ch_cary <- function(grid, path = ""){
  file <- "/projectnb/dietzelab/fosterj/Data/Cary_mouse.csv"
  
  dat <- read.csv(paste(path, file, sep = ""),
                  na.strings=c(c("", " ", "      "),"NA")) # read in data
  
  smam <- dat[,c("Grid",
                 "Full.Date.1",
                 "Full.Date.2",
                 "Day.1",
                 "Day.2",
                 "Tag..",
                 "Fate")]
  
  alive <- c(1, 2)                                    # codes for alive individuals
  #smam <- smam %>% filter(Fate %in% alive)            # extract only live individuals
  smam <- subset(smam, Fate == 1 | Fate == 2)
  
  smam[4:5] <- as.integer(!is.na(smam[4:5]))          # converts trap histories to 1 or 0
  
  smam$Full.Date.1 <- as.Date(smam$Full.Date.1)           # convert factors to dates
  smam$Full.Date.2 <- as.Date(smam$Full.Date.2)           # convert factors to dates
  
  m <- subset(smam, Grid == grid)
  m <- m[,c("Tag..", "Day.1", "Day.2", "Full.Date.1", "Full.Date.2")]
  m <- subset(m, !is.na(Tag..))
  
  m$Tag.. <- as.character(m$Tag..)
  m.ls <- split(m, m$Full.Date.1)                     # split on first capture date for sampling occasion
  for(i in 1:length(m.ls)){
    m.ls[[i]] <- m.ls[[i]][c(-4, -5)]
  }
  
  day.1 <- as.character(unique(m$Full.Date.1))        # 1st capture date of sampling occasion
  day.2 <- as.character(unique(m$Full.Date.2))        # 2nd capture date of sampling occasion
  days <- c(rbind(day.1, day.2))                      # vector of unique trapping days (for colnames)
  
  ch.base <- merge(m.ls[[1]], m.ls[[2]], by = "Tag..", all = TRUE)
  l.m <- length(m.ls)*1
  for(i in 3:l.m){                                    # loop through the rest
    g <- as.data.frame(m.ls[[i]])
    ch.base <- merge(ch.base, g, by = "Tag..", all = TRUE)
  }
  ch.base <- as.matrix(ch.base[,-1])                            # convert all NAs to 0
  for(i in 1:nrow(ch.base)){
    for(t in 1:ncol(ch.base)){
      if(is.na(ch.base[i,t])){ch.base[i,t]<-0}
      as.numeric(ch.base[i,t])
    }
  }
  colnames(ch.base) <- days
  return(ch.base)
}