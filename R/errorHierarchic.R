##' Later
##' @param nlnet The nlnet
##' @param trainIn training data 
##' @param trainOut  fitted data
##' @return error
##' @author Henning Redestig, Matthias Scholz
errorHierarchic <- function(nlnet, trainIn, trainOut) {

  weights <- nlnet@weights$current()
  
  if(nlnet@inverse) {
    numElements <- nlnet@net[1] * dim(trainOut)[2]
    trainIn <- array(weights[1:numElements], dim=c(nlnet@net[1], dim(trainOut)[2]))
    wTrainIn <- weights[1:numElements, drop=FALSE]
    weights <- weights[(numElements + 1):length(weights),,drop=FALSE]
  }

  netDim <- dim(nlnet@net)
  trainDim <- dim(trainOut)
  weightMats <- vector2matrices(weights, nlnet@net)

  hierarchicIdx <- nlnet@hierarchic$idx[,nlnet@hierarchic$var != 0, drop=FALSE]
  hierarchicVar <- nlnet@hierarchic$var[,colSums(nlnet@hierarchic$var) != 0, drop=FALSE]

  subnetNum <- length(hierarchicVar)

  out <- array(0, dim=c(trainDim[1], trainDim[2], subnetNum))

  sBias <- array(1, dim=c(1, trainDim[2]))

  sExtract <- eval(parse(text=paste(nlnet@fct[1], "(trainIn)")))
  
  Eitemize <- NULL
  if(nlnet@hierarchic$layer > 1) {      #this should not be executed at all if sequence is 1:0
    for(layer in 1:(nlnet@hierarchic$layer - 1)) {
      sExtract <- rbind(sBias, sExtract)
      sExtract <- eval(parse(text=paste(nlnet@fct[layer + 1], "(weightMats[[layer]] %*% sExtract)")))
    }
  }

  for(subnet in 1:subnetNum) {
    sRecon <- sExtract
    sRecon[hierarchicIdx[,subnet]==0,] <- 0

    for(layer in nlnet@hierarchic$layer:(netDim[2] - 1)) {
      sRecon <- rbind(sBias, sRecon)
      sRecon <- eval(parse(text=paste(nlnet@fct[layer+1], "(weightMats[[layer]] %*% sRecon)")))
    }
    out[,,subnet] <- sRecon

    ## error function
    eTmp <- (sRecon - trainOut)^2
    eTmp[is.na(eTmp)] <- 0
    
    Eitemize[subnet] <- sum(eTmp) * 0.5

    if(!is.null(nlnet@dataDist)) 
      Eitemize[subnet] <- 0.5 * sum(nlnet@dataDist * eTmp)
    else
      Eitemize[subnet] <- 0.5 * sum(eTmp)
  }

  error <- tcrossprod(hierarchicVar, rbind(Eitemize))

  if(!is.null(nlnet@weightDecay))
    error <- error + nlnet@weightDecay * 0.5 * sum(weights^2)

  ## smooth (0.01) weight decay also for input values
  if(nlnet@inverse)
    error <- error + 0.01 * nlnet@weightDecay * 0.5 * sum(wTrainIn^2)
  
  return(list(error=error, out=out))
}
  

