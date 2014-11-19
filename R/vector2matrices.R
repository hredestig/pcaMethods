##' Tranform the vectors of weights to matrix structure
##' @param object  an nlpcaNet
##' @return weights in matrix structure
##' @author Henning Redestig
##' @aliases vector2matrices,nlpcaNet-method
setMethod("vector2matrices", "nlpcaNet", function(object) {
  netDim <- dim(object@net)
  posBegin <- 1
  posEnd <- 0
  result <- list()
  for(i in 1:(netDim[2] - 1)) {
    wSize <- c(object@net[i + 1], object@net[i] + 1)
    posEnd <- posEnd + prod(wSize)
    result[[i]] <-
      matrix(object@weights$current()[posBegin:posEnd], wSize[1], wSize[2])
    posBegin <- posEnd + 1
  }
  
  if(posEnd < length(object@weights$current()))
    stop("weight vector has too many elements\n")
  result
})

##' Tranform the vectors of weights to matrix structure
##' @param object an nlpcaNet
##' @param net the neural network
##' @return weights in matrix structure
##' @author Henning Redestig
##' @aliases vector2matrices,matrix-method
setMethod("vector2matrices", "matrix", function(object, net) {
  netDim <- dim(net)
  posBegin <- 1
  posEnd <- 0
  result <- list()
  for(i in 1:(netDim[2] - 1)) {
    wSize <- c(net[i + 1], net[i] + 1)
    posEnd <- posEnd + prod(wSize)
    result[[i]] <- matrix(object[posBegin:posEnd], wSize[1], wSize[2])
    posBegin <- posEnd + 1
  }
  if(posEnd < length(object))
    stop("weight vector has too many elements\n")
  result
})
