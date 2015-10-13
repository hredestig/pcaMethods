##' Conjugate gradient optimization
##' @param nlnet The nlnet
##' @param trainIn Training data
##' @param trainOut fitted data
##' @param verbose logical, print messages
##' @return ...
##' @author Henning Redestig, Matthias Scholz
optiAlgCgd <- function(nlnet,
                       trainIn,
                       trainOut,
                       verbose=FALSE) {
  
  tmpnet <- forkNlpcaNet(nlnet)
  derr <- tmpnet@gradient(tmpnet, trainIn, trainOut)
  dw <- derr$dwTotal
  e <- derr$Etotal

  dv <- -dw

  if(tmpnet@featureSorting)
    eSortLast <- e

  eHist <- rep(0, tmpnet@maxIter)
  ttLast <- rep(0.0001, 6)

  for(i in 1:tmpnet@maxIter) {
    if(verbose) {
      if(i %% 10 == 0)
        cat("*")
      if(i %% 100 == 0)
        cat(" [", i, "]\n")
    }

    eHist[i] <- e
    eLast <- e
    # line search in direction dv (downhill)
    ttGuess <- max(min(ttLast), 0.00001)
    
    linSe <- lineSearch(tmpnet, dv, e, ttGuess, trainIn, trainOut, verbose)
    tmpnet@weights$set(cbind(linSe$wBest))
    e <- linSe$eBest
    tt <- linSe$ttBest
    
    ttLast <- c(ttLast[2:length(ttLast)], tt) #shift and add new tt

    gradRes <- tmpnet@gradient(tmpnet, trainIn, trainOut)
    dwNew <- gradRes$dwTotal
    e <- gradRes$Etotal

    ## define new search direction dv (conjugate direction)
    
    ## b1=dw_new'*dw_new; # Fletcher-Reeves
    b1 <- crossprod(dwNew, (dwNew - dw))#Polak-Ribiere
    b2 <- crossprod(dw)
    beta <- b1 / b2
    
    dv <- -dwNew + dv %*% beta
    dw <- dwNew

    if(e > eLast) {
      dv <- -dwNew
      if(verbose)
        cat("!", sep="")
    }

    if(is.na(e))
      stop("Square error is NA (critical) - accuracy in line-search might be too small")

    if(tmpnet@featureSorting) 
      if(e / eSortLast < 0.90 || i == tmpnet@maxIter ||
         i == tmpnet@maxIter - 1 || i == tmpnet@maxIter - 2) {
        eSortLast <- e
        if(verbose)
          cat("<", i, ">", sep="")
        ## somewhat secret method, sortFeatures calls
        ## nlnet@weights$set(x) so the weights are updated here
        ## 'behind the scenes'
        sortFeatures(tmpnet, trainIn, trainOut)
      }

  }
  tmpnet
}
