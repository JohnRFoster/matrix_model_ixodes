##' Continues to run MCMC until convergence and large enough effective sample size
##'
##' @param j.model The jags model
##' @param variableNames A vector of the variable names that need to converge
##' @param maxIter The maximum number of iterations to run
##' @param baseNum The number of initial iterations to run
##' @param iterSize The number of iterations per each run (between when GBR values are checked)
##' @param maxGBR The maximum allowable GBR value after the baseNum of iterations are run
##' @param ID An optional identifier to print if it doesn't converge
##' @param sampleCutoff The minimum desired effective sample size for testing (default is 5000)
##' @import rjags
##' @import runjags
##' @import coda
##' @export

library(runjags)

runMCMC_Model <- function(j.model, variableNames,
                          tick = FALSE,
                          mice = FALSE,
                          maxIter = 1000000000,
                          baseNum = 100000,
                          iterSize = 50000,
                          maxGBR = 13,
                          ID = "",
                          sampleCutoff = 3000,
                          ...){
  
  first.call <- Sys.time()
  cat("\nInitial coda.samples call \n")
  var.out   <- coda.samples (model = j.model,
                             variable.names = variableNames,
                             n.iter = baseNum,
                             thin = ceiling(baseNum / 5000),
                             ...)
  
  time.elapsed <- Sys.time() - first.call 
  cat(baseNum, "iterations after:\n")
  print(time.elapsed)
  
  numb <- baseNum
  continue <- TRUE
  GBR.bad <- TRUE
  burnin <- 0
  while(continue & numb<maxIter){
    loop.start <- Sys.time()
    new.out   <- coda.samples(model = j.model,
                              variable.names = variableNames,
                              n.iter = iterSize,
                              thin = ceiling(iterSize / 5000))
    
    var.out <- combine.mcmc(mcmc.objects=list(var.out,new.out),collapse.chains = FALSE)
    params <- var.out
    
    time.elapsed <- Sys.time() - loop.start 
    cat(iterSize, "current iterations in:\n")
    print(time.elapsed)
    
    time.elapsed <- Sys.time() - first.call 
    cat(numb, "iterations complete for a total run time of:\n")
    print(time.elapsed)
    
    if(tick){
      params <- list()
      mfit <- as.matrix(var.out, chains = TRUE)
      pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
      m.cols <- grep("m[", colnames(mfit), fixed = TRUE)
      chain.col <- which(colnames(mfit) == "CHAIN")
      params <- ecoforecastR::mat2mcmc.list(mfit[, -c(m.cols, pred.cols)])  
    }
    
    continue <- FALSE
    if(GBR.bad){
      GBR.vals <- gelman.diag(params, multivariate = FALSE)
      GBR.bad <- FALSE
      for(i in 1:nrow(GBR.vals$psrf)){
        for(j in 1:ncol(GBR.vals$psrf)){
          if(GBR.vals$psrf[i,j]>maxGBR){
            print(GBR.vals)
            cat("\n GBR values too high:", ID)
            return(FALSE)
          }
          if(GBR.vals$psrf[i,j]>1.1){
            continue <-  TRUE
            GBR.bad <- TRUE
          }
        }
      }
      print(GBR.vals)
    }
    if(!continue){
      if(burnin==0){
        GBR <- gelman.plot(params)
        # burnin <- GBR$last.iter[tail(which((GBR$shrink[,,2]>1.05)),1)+1]
        burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.05,1,any)),1)+1]
        if(length(burnin) == 0) burnin = 1
      }
      message("Burnin: ", burnin)
      var.burn <- window(var.out,start=burnin)
      #var.burn <- var.out
      effsize <- effectiveSize(params)
      for(i in 1:length(effsize)){
        if(effsize[i]<sampleCutoff){
          continue = TRUE
        }
      }
      cat("\nEffective size:", effsize, "\n")
      # save(var.burn)
    }
    numb <- numb+iterSize
  }
  if(continue==TRUE){
    cat("\n Model Did not Converge")
    var.burn <- FALSE
  }
  return(var.burn)
}

