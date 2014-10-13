##' Later
##' @param nlnet the nlnet
##' @param trainIn training data
##' @param trainOut  fitted data
##' @return derror
##' @author Henning Redestig, Matthias Scholz
derrorHierarchic <- function(nlnet, trainIn, trainOut) {

  weights <- nlnet@weights$current()
  netDim <- dim(nlnet@net) 
  if(nlnet@inverse) {
    numElements <- nlnet@net[1] * dim(trainOut)[2]
    trainIn <- matrix(weights[1:numElements], nrow=nlnet@net[1], ncol=dim(trainOut)[2])
    wTrainIn <- weights[1:numElements,drop=FALSE]
    weights <- weights[(numElements + 1):length(weights), ,drop=FALSE]
  }

  weightMats <- vector2matrices(weights, nlnet@net)
  trainDim <- dim(trainIn)
  subnetNum <- length(nlnet@hierarchic$var)

  ## ******************************

  Epattern <- array(0, dim=c(dim(trainOut), subnetNum))
  nOut <- array(0, dim=c(sum(nlnet@net), trainDim[2], subnetNum))

  for(subnet in 1:subnetNum)
    nOut[1:trainDim[1],,subnet] <- eval(parse(text=paste(nlnet@fct[1], "(trainIn)")))

  if(nlnet@inverse)
    for(subnet in 1:subnetNum)
      nOut[nlnet@hierarchic$idx[,subnet]==0,,subnet] <- 0

  ## forward propagation
  for(subnet in 1:subnetNum) {
    if(nlnet@hierarchic$var[subnet] != 0) {
      sBias <- array(1, dim=c(1, trainDim[2]))
      for(i in 1:(netDim[2] - 1)) {
        if(i == 1)
          nBegin <- 1
        else
          nBegin <- sum(nlnet@net[1:(i-1)])+1
        sIn <- rbind(sBias, nOut[nBegin:sum(nlnet@net[1:i]),,subnet])

        sOut <- eval(parse(text=paste(nlnet@fct[i+1], "(weightMats[[i]] %*% sIn)")))

        if(i == (nlnet@hierarchic$layer - 1)) 
          sOut[nlnet@hierarchic$idx[,subnet]==0,] <- 0

        nOut[(sum(nlnet@net[1:i])+1):sum(nlnet@net[1:(i+1)]),,subnet] <- sOut
      }
      output <- nOut[(sum(nlnet@net[1:(length(nlnet@net)-1)])+1):dim(nOut)[1],,subnet]
      Epattern[,,subnet] <- output - trainOut
    }
  }

  ## error function
  Epattern <- Epattern^2
  Epattern[is.na(Epattern)] <- 0        #set the missing values to zero

  if(!is.null(nlnet@dataDist))
    for(subnet in 1:subnetNum)
      Epattern[,,subnet] <- Epattern[,,subnet] * nlnet@dataDist

  Eitemize <- apply(Epattern, 3, sum) * 0.5
  Etotal <- sum(nlnet@hierarchic$var * Eitemize)

  if(!is.null(nlnet@weightDecay))
    Etotal <- Etotal + nlnet@weightDecay * 0.5 * sum(weights^2)

  if(nlnet@inverse)
    Etotal <- Etotal + 0.01 * nlnet@weightDecay * 0.5 * sum(wTrainIn^2)

  ## back propagation
  nError <- array(0, dim=c(sum(nlnet@net), trainDim[2], subnet))
  dWeight <- vector(length=netDim[2] - 1, mode="list")
  wBp <- vector(length=netDim[2] - 1, mode="list")

  ## wBp is weights for back propagation
  for(u in 1:(netDim[2] - 1)) 
    wBp[[u]] <- weightMats[[u]][,2:(nlnet@net[u] + 1)] # cats the weights which belong to bias

  dw <- array(0, dim=c(length(weights), subnet))
  for(subnet in 1:subnetNum) {
    if(nlnet@hierarchic$var[subnet] != 0) {
      ## last layer
      sTmp <- nOut[(dim(nOut)[1]-nlnet@net[length(nlnet@net)]+1):dim(nOut)[1],,subnet]
      if(nlnet@fct[length(nlnet@fct)] == "tanh")
        eTmp <- (1 - sTmp^2) * (sTmp - trainOut) #prev trainOut - sTmp (fixed to get rid of sign change)
      else if(nlnet@fct[length(nlnet@fct)] == "linr")
        eTmp <- sTmp - trainOut #prev trainOut - sTmp (fixed to get rid of sign change)
      eTmp[is.na(eTmp)] <- 0
      if(!is.null(nlnet@dataDist))
        eTmp <- eTmp * nlnet@dataDist
      
      nError[(dim(nError)[1]-nlnet@net[length(nlnet@net)]+1):dim(nError)[1],,subnet] <- eTmp

      ## all other layers
      for(n in 1:(netDim[2] - 1)){
        i <- netDim[2]-n
        
        ## the if clause is to avoid 1:0 difference in R
        ## Matlab (1:0 => Empty matrix), R (1:0 => [1,0])
        if(i > 1)
          sTmp <- nOut[(sum(nlnet@net[1:(i-1)])+1):sum(nlnet@net[1:i]),,subnet]
        else
          sTmp <- nOut[1:sum(nlnet@net[1:i]),,subnet]

        if(i==(nlnet@hierarchic$layer-1))
          eTmp[nlnet@hierarchic$idx[,subnet]==0,] <- 0

        dWeight[[i]] <- tcrossprod(eTmp, rbind(sBias, sTmp)) #gradient

        if (nlnet@fct[i] == "tanh")
          eTmp <- (1 - sTmp^2) * crossprod(wBp[[i]],eTmp)
        else if (nlnet@fct[i] == "linr")
          eTmp <- crossprod(wBp[[i]], eTmp)

        ## the if clause is to avoid 1:0 difference in R
        if(i > 1) 
          nError[(sum(nlnet@net[1:(i - 1)]) + 1):sum(nlnet@net[1:i]), ,subnet] <- eTmp
        else 
          nError[1:sum(nlnet@net[1:i]), ,subnet] <- eTmp
      }
      dw[,subnet] <- unlist(dWeight)    #fixed sign change
      
    }
  }

  if(nlnet@inverse) {
    dw <- rbind(array(0, dim=c(numElements, subnetNum)), dw)
    for(subnet in 1:subnetNum) {
      eTmp <- array(nError[1:nlnet@net[1],,subnet], dim=c(nlnet@net[1], dim(nError)[2]))
      
      eTmp[nlnet@hierarchic$idx[,subnet] == 0,] <- 0

      dim(eTmp) <- NULL                 #a bit unsure if this is correct but seems to work
      dw[1:numElements,subnet] <- unlist(eTmp) #fixed sign change
    }
    ## weights <- rbind(cbind(rep(0, numElements)), cbind(weights)) #old: only weight decay for real weights
    weights <- rbind(cbind(0.01 * wTrainIn), cbind(weights)) #new 
  }

  dwTotal <- array(0, dim=dim(weights))
  for (subnet in 1:subnetNum)  {
    dwTotal <- dwTotal + nlnet@hierarchic$var[subnet] * dw[, subnet]
  }

  if(!is.null(nlnet@weightDecay))
    dwTotal <- dwTotal + nlnet@weightDecay * weights
  return(list(dwTotal=dwTotal, Etotal=Etotal, nError=nError, nOut=nOut))
}

