##---------------------------------------------------------------------------------##
##                    Known states for capture history matrix                      ##
##---------------------------------------------------------------------------------##

#' Known states function
#'
#' This function fills in the capture history matrix where individuals are known but not observed
#' @param ch Capture history matrix
#' @export
#' @examples known_states(ch)

known_states <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] != 0))
    n2 <- max(which(ch[i,] != 0))
    if(n2 > n1){state[i, n1:n2] <- 1}
    #state[i, n1:n2] <- 1
  }
  return(state)
}