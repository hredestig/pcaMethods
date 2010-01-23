##' Sort the features of NLPCA object
##' @param nlnet The nlnet
##' @param trainIn Training data in
##' @param trainOut  Training data after it passed through the net
##' @return ...
##' @author Henning Redestig
sortFeatures <- function(nlnet, trainIn, trainOut) {
  weightsAll <- nlnet@weights$current()
  weights <- weightsAll
  if(nlnet@inverse) {
    numElements <- nlnet@net[1] * dim(trainOut)[2]
    trainIn <- array(unlist(weightsAll), dim=c(nlnet@net[1], dim(trainOut)[2]))
    weights <- weightsAll[(numElements + 1):length(weightsAll),,drop=FALSE]
  }
  
  netDim <- dim(nlnet@net)
  trainDim <- dim(trainIn)
  bneckNum <- nlnet@net[nlnet@componentLayer]
  weightMats <- vector2matrices(weights, nlnet@net)

  bneckNum <- nlnet@net[nlnet@componentLayer]

  ## ******************************
  nOut <- array(0, dim=c(sum(nlnet@net), trainDim[2], 2))
  for(subnet in 1:2)
    nOut[1:trainDim[1],,subnet] <- trainIn
  
  ## forward propagation
  for(n in 0:(bneckNum - 2)) {
    E <- c(0,0)
    for(choice in 1:2) {
      sBias <- rep(1, trainDim[2])
      for(i in 1:(netDim[2] - 1)) {
        if(i == 1)
          nBegin <- 1
        else
          nBegin <- sum(nlnet@net[1:(i - 1)]) + 1
        sIn <- rbind(sBias, nOut[nBegin:sum(nlnet@net[1:i]),, choice])
        sOut <- eval(parse(text=paste(nlnet@fkt[i], "(weightMats[[i]] %*% sIn)")))
        
        if(i == nlnet@componentLayer - 1) {
          idx <- rep(0, bneckNum)
          idx[1:(n + choice)] <- 1
          
          if(choice == 2)
            idx[n+choice-1] <- 0
          sOut[idx == 0,] <- 0
        }
        
        nOut[(sum(nlnet@net[1:i]) + 1):sum(nlnet@net[1:(i+1)]),,choice] <- sOut
      }
      output <-
        nOut[(sum(nlnet@net[1:(dim(nlnet@net)[2]-1)])+1):dim(nOut)[1], ,choice]
      
      Epattern <- (output - trainOut)^2
      Epattern[is.na(Epattern)] <- 0
      if(!is.null(nlnet@dataDist))
        Epattern <- Epattern * nlnet@dataDist
      E <- mean(Epattern)
      E[choice] <- E
    }

    if(E[1]>E[2]) {                     #change features
      changeIdx <- 1:bneckNum
      changeIdx[(n+1):(n+2)] <- c(n+2, n+1)
      weightMats[[nlnet@componentLayer - 1]] <-
        weightMats[[nlnet@componentLayer - 1]][changeIdx,]
      weightMats[[nlnet@componentLayer]] <-
        weightMats[[nlnet@componentLayer]][,c(1,changeIdx+1)]
      switching <- c(n+1, n+2)
      nlnet@fCount <- as.integer(nlnet@fCount + 1)
    }
  }
  weights <- cbind(unlist(weightMats))

  if(nlnet@inverse)
    nlnet@weights$set(rbind(matrix(trainIn, nrow=numElements, ncol=1), weights))
}
