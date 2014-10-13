##' Neural network based non-linear PCA
##'
##' Artificial Neural Network (MLP) for performing non-linear
##' PCA. Non-linear PCA is conceptually similar to classical PCA but
##' theoretically quite different. Instead of simply decomposing our
##' matrix (X) to scores (T) loadings (P) and an error (E) we train a
##' neural network (our loadings) to find a curve through the
##' multidimensional space of X that describes a much variance as
##' possible. Classical ways of interpreting PCA results are thus not
##' applicable to NLPCA since the loadings are hidden in the network.
##' However, the scores of components that lead to low
##' cross-validation errors can still be interpreted via the score
##' plot.  Unfortunately this method depend on slow iterations which
##' currently are implemented in R only making this method extremely
##' slow. Furthermore, the algorithm does not by itself decide when it
##' has converged but simply does 'maxSteps' iterations.
##' @title Non-linear PCA
##' @param Matrix \code{matrix} --- Preprocessed data with the
##' variables in columns and observations in rows. The data may
##' contain missing values, denoted as \code{NA}
##' @param nPcs \code{numeric} -- Number of components to
##' estimate. The preciseness of the missing value estimation depends
##' on thenumber of components, which should resemble the internal
##' structure of the data.
##' @param maxSteps \code{numeric} -- Number of estimation
##' steps. Default is based on a generous rule of thumb.
##' @param unitsPerLayer The network units, example: c(2,4,6) for two
##' input units 2feature units (principal components), one hidden
##' layer fornon-linearity and three output units (original amount
##' ofvariables).
##' @param functionsPerLayer The function to apply at each layer
##' eg. c("linr", "tanh", "linr") 
##' @param weightDecay Value between 0 and 1.
##' @param weights Starting weights for the network. Defaults to
##' uniform random values but can be set specifically to make
##' algorithm deterministic.
##' @param verbose \code{boolean} -- nlpca prints the number of steps
##' and warning messages if set to TRUE. Default is interactive().
##' @param ...  Reserved for future use. Not passed on anywhere.
##' @return Standard PCA result object used by all PCA-basedmethods of
##' this package. Contains scores, loadings, data meanand more. See
##' \code{\link{pcaRes}} for details.
##' @author Based on a matlab script by Matthias Scholz and ported to
##' R by Henning Redestig
##' @references Matthias Scholz, Fatma Kaplan, Charles L Guy, Joachim
##' Kopkaand Joachim Selbig. Non-linear PCA: a missing
##' data approach. \emph{Bioinformatics, 21(20):3887-3895, Oct 2005}
##' @examples
##' ## Data set with three variables where data points constitute a helix
##' data(helix)
##' helixNA <- helix
##' ## not a single complete observation
##' helixNA <- t(apply(helix, 1, function(x) { x[sample(1:3, 1)] <- NA; x}))
##' ## 50 steps is not enough, for good estimation use 1000
##' helixNlPca <- pca(helixNA, nPcs=1, method="nlpca", maxSteps=50)
##' fittedData <- fitted(helixNlPca, helixNA)
##' plot(fittedData[which(is.na(helixNA))], helix[which(is.na(helixNA))])
##' ## compared to solution by Nipals PCA which cannot extract non-linear patterns
##' helixNipPca <- pca(helixNA, nPcs=2)
##' fittedData <- fitted(helixNipPca)
##' plot(fittedData[which(is.na(helixNA))], helix[which(is.na(helixNA))])
##' @export
nlpca <- function(Matrix, nPcs=2, 
                  maxSteps=2 * prod(dim(Matrix)),
                  unitsPerLayer=NULL, functionsPerLayer=NULL,
                  weightDecay=0.001, weights=NULL,
                  verbose=interactive(),...) {

  ## do some basic checks
  object <- Matrix
  trainIn <- NULL
  trainOut <- t(object)
  stds <- apply(trainOut, 2, sd, na.rm=TRUE)
  scalingFactor <- 0.1 / max(stds)
  trainOut <- trainOut * scalingFactor

  ## now setup the initial nlpcaNet object
  numNaN <- sum(is.na(object))
  
  ## always inverse in this version, bottleneck is not fully implemented 
  inverse <- TRUE

  ## DATADIST (nlnet@dataDist) is given by weightOut
  dataDist <- apply(!is.na(trainOut), 2, as.integer) #0 for NA, 1 for everything else
  if(!inverse)
    dataDist <- NULL

  ## setup the network architecture
  if(is.null(unitsPerLayer)) {
    ld <- dim(trainOut)[1]
    lh <- nPcs
    if(nPcs < 10)
      lh <- 2 + 2 * nPcs
    unitsPerLayer <- c(ld, lh, nPcs, lh, ld)
    if(inverse)
      unitsPerLayer <- c(nPcs, lh, ld)
  }
  featureLayer <- ceiling(length(unitsPerLayer) / 2)
  if(inverse)
    featureLayer <- 1

  if(is.null(functionsPerLayer)) {
    functionsPerLayer <- rep("tanh", length(unitsPerLayer))
    functionsPerLayer[1] <- "linr"
    functionsPerLayer[featureLayer] <- "linr"
    functionsPerLayer[length(unitsPerLayer)] <- "linr"
  }
  hierarchic <- list(layer=featureLayer,
                     var=rbind(c(rep(1, nPcs), 0.01)),
                     idx=getHierarchicIdx(unitsPerLayer[featureLayer]))
  
  ## set up the weights
  wNum <- sum(sapply(2:length(unitsPerLayer),
                     function(i) (1 + unitsPerLayer[i - 1]) * unitsPerLayer[i]))
  if(!is.null(weights) && length(weights) != wNum) {
    warning("Weight vector not expected length (", wNum, "), using random weights", sep="")
    weights <- NULL
  }
  if(is.null(weights)) 
    weights <- cbind(0.2 * (runif(wNum, 0, 1) - 0.1))
  
  if(inverse) {
    numPattern <- dim(trainOut)[2]
    tmpTrainIn <- cbind(rnorm(unitsPerLayer[1] * numPattern,0,1) * 0.1)
    weights <- rbind(tmpTrainIn, weights)
  }

  if(nPcs == 1)
      featureSorting <- FALSE
  if(nPcs > 1)
    featureSorting <- TRUE

  nlnet <- new("nlpcaNet")
  nlnet@net <- rbind(unitsPerLayer)
  nlnet@hierarchic <- hierarchic
  nlnet@fct <- functionsPerLayer
  nlnet@fkt <- functionsPerLayer[2:length(functionsPerLayer)]
  nlnet@weightDecay <- weightDecay
  nlnet@featureSorting <- featureSorting
  nlnet@dataDist <- dataDist
  nlnet@inverse <- inverse
  nlnet@fCount <- as.integer(0)
  nlnet@componentLayer <- as.integer(featureLayer)
  nlnet@error <- errorHierarchic
  nlnet@gradient <- derrorHierarchic
  nlnet@maxIter <- as.integer(maxSteps)
  nlnet@weights <- weightsAccount(weights)
  nlnet@scalingFactor <- scalingFactor
  ## ******************************

  if(verbose)
    cat("Training network with", nlnet@maxIter, "iterations...\n!:\tSquare error is NA -- accuracy in line-search might be too small\n<n>:\tComponents were sorted at iteration n\n^:\tToo many iterations while expanding\n")
  newnet <- optiAlgCgd(nlnet, trainIn, trainOut, verbose)
  if(verbose)
    cat("\nDone\n")

  if(inverse) {
    nObs <- unitsPerLayer[1] * dim(trainOut)[2]
    we <- newnet@weights$current()
    scores <- t(matrix(we[1:nObs], nrow=unitsPerLayer[1], dim(trainOut)[2]))
    newnet@weights$set(we[(nObs + 1):length(we),,drop=FALSE])
  }
  ## for further applications newnet must not be inverse anymore
  newnet@inverse <- FALSE               

  res <- new("pcaRes")
  res@scores <- scores
  res@loadings <- matrix()
  res@network <- newnet
  res@method <- "nlpca"

  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for(i in 1:nPcs) 
    R2cum[i] <-
      1 - sum((Matrix - fitted(res, Matrix, nPcs=i))^2, na.rm=TRUE) / TSS

  res@R2cum <- R2cum
  res
}

##' Index in hiearchy
##' @param hierarchicNum  A number
##' @return ...
##' @author Henning Redestig, Matthias Scholz
getHierarchicIdx <- function(hierarchicNum) {
  res <- matrix(1, ncol=hierarchicNum, nrow=hierarchicNum)
  res[lower.tri(res)] <- 0
  cbind(res, c(0, rep(1, hierarchicNum - 1)))
}

