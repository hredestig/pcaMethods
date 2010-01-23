##' Complete copy of nlpca net object
##' @param nlnet a nlnet
##' @return A copy of the input nlnet
##' @author Henning Redestig
forkNlpcaNet <- function(nlnet) {
  res <- new("nlpcaNet")
  res@net <- nlnet@net
  res@hierarchic <- nlnet@hierarchic
  res@fct <- nlnet@fct
  res@fkt <- nlnet@fkt
  res@weightDecay <- nlnet@weightDecay
  res@featureSorting <- nlnet@featureSorting
  res@dataDist <- nlnet@dataDist
  res@inverse <- nlnet@inverse
  res@fCount <- nlnet@fCount
  res@componentLayer <- nlnet@componentLayer
  res@error <- nlnet@error
  res@gradient <- nlnet@gradient
  res@weights <- weightsAccount(nlnet@weights$current())
  res@maxIter <- nlnet@maxIter
  res@scalingFactor <- nlnet@scalingFactor
  res
}
