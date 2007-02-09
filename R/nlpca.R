nlpca <- function(Matrix, nPcs=2, center=TRUE, completeObs=TRUE, maxSteps=2 * prod(dim(Matrix)),
                  unitsPerLayer=NULL, functionsPerLayer=NULL, weightDecay=0.001, weights=NULL,
                  verbose=interactive(),...) {

  ##<..Beg Rdocu..>
  ## ~name~
  ##   nlpca
  ## ~title~
  ##   NLPCA
  ## ~description~
  ##   Neural network based non-linear PCA
  ## ~usage~
  ##   nlpca(Matrix, nPcs=2, center=TRUE, completeObs=TRUE, maxSteps=2*prod(dim(Matrix)), unitsPerLayer=NULL, functionsPerLayer=NULL, weightDecay=0.001, weights=NULL, verbose=interactive(), ...)
  ## ~arguments~
  ##   ~-Matrix~
  ##     \code{matrix} --- Data containing the variables in columns
  ##     and observations in rows. The data may contain missing
  ##     values, denoted as \code{NA}
  ##   ~-nPcs~
  ##     \code{numeric} -- Number of components to estimate. The
  ##     preciseness of the missing value estimation depends on the
  ##     number of components, which should resemble the internal
  ##     structure of the data.
  ##   ~-center~
  ##     \code{boolean} Mean center the data if TRUE
  ##   ~-completeObs~
  ##     \code{boolean} Return the complete observations if TRUE. This
  ##     is the original data with NA values filled with the estimated
  ##     values.
  ##   ~-maxSteps~
  ##     \code{numeric} -- Number of estimation steps. Default is based on a generous heuristic.
  ##   ~-unitsPerLayer~
  ##     The network units, example: c(2,4,6) for two input units 2
  ##     feature units (principal components), one hidden layer for
  ##     non-linearity and three output units (original amount of
  ##     variables).
  ##   ~-functionsPerLayer~
  ##     The function to apply at each layer eg. c("linr", "tanh", "linr") 
  ##   ~-weightDecay~
  ##     Value between 0 and 1.
  ##   ~-weights~
  ##     Starting weights for the network. Defaults to uniform random
  ##     values but can be set specifically to make algorithm
  ##     deterministic.
  ##   ~-verbose~
  ##     \code{boolean} -- BPCA prints the number of steps and the
  ##     increase in precision if set to TRUE. Default is
  ##     interactive().
  ##   ~-...~
  ##     Reserved for future use. Not passed on anywhere.
  ## ~details~
  ##   Artificial Neural Network (MLP) for performing non-linear PCA. 
  ## ~value~
  ##   \item{pcaRes}{Standard PCA result object used by all PCA-based
  ##   methods of this package. Contains scores, loadings, data mean
  ##   and more. See \code{\link{pcaRes}} for details.
  ## ~references~
  ##    Matthias Scholz, Fatma Kaplan, Charles L Guy, Joachim Kopka
  ##    and Joachim Selbig. Non-linear PCA: a missing data
  ##    approach. \emph{Bioinformatics, 21(20):3887-3895, Oct 2005}
  ## ~examples~
  ##   data(helix)
  ##   helixNA <- helix
  ##   
  ##   helixNA[sample(1:3, 1000, replace=TRUE)] <- NA # not a single complete observation
  ##   
  ##   helixPca <- pca(helixNA, nPcs=1, method="nlpca", maxSteps=500)
  ##   
  ##   plot(helixPca@completeObs[which(is.na(helixNA))], helix[which(is.na(helixNA))])
  ## ~keywords~
  ##   multivariate
  ## ~author~
  ##   Based on a matlab script by Matthias Scholz 
  ##   <matthias.scholz[at]uni-greifswald.de> and ported to R by Henning
  ##   Redestig <redestig[at]mpimp-golm.mpg.de>
  ##>..End Rdocu..<
  
  ## debug {
  ##   Matrix <- t(dataNA)
  ##   nPcs=2; center=TRUE; completeObs=TRUE; maxSteps=2000;
  ##   functionsPerLayer=NULL; weightDecay=0.001; weights=NULL;
  ##   verbose=interactive(); unitsPerLayer <- NULL; 
  ## }

  ## do some basic checks
  if (!checkData(Matrix, verbose=verbose))
    stop("Invalid data format, use checkData(Matrix, verbose = TRUE) for details.\n")
  if (nPcs > ncol(Matrix))
    stop("more components than matrix columns selected, exiting")

  if (center) {
    object <- scale(Matrix, center = TRUE, scale = FALSE)
    means <- attr(object, "scaled:center")
  } else
  object <- Matrix
  trainIn <- NULL
  trainOut <- t(object)
  stds <- apply(trainOut, 2, sd, na.rm=TRUE)
  scalingFactor <- 0.1 / max(stds)
  trainOut <- trainOut * scalingFactor

  ## ******************************
  ## now setup the initial nlpcaNet object
  ## ******************************
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

  if (completeObs) {
    Ye <- errorHierarchic(newnet, trainIn, trainOut)$out[,,nPcs]
    Ye <- Ye / scalingFactor
    if (center) 
        Ye <- Ye + means
    cObs <- t(object)
    cObs[is.na(cObs)] <- Ye[is.na(cObs)]
  }

  if(inverse) {
    nObs <- unitsPerLayer[1] * dim(trainOut)[2]
    we <- newnet@weights$current()
    scores <- t(matrix(we[1:nObs], nrow=unitsPerLayer[1], dim(trainOut)[2]))
    newnet@weights$set(we[(nObs + 1):length(we),,drop=FALSE])
  }
  newnet@inverse <- FALSE               #for further applications newnet must not be inverse anymore

  
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- t(cObs)
  r@scores <- scores
  r@network <- newnet
  r@sDev <- apply(scores, 2, sd)
  r@nObs <- nObs
  r@nVar <- ncol(object)
  r@nPcs <- nPcs
  r@centered <- center
  r@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
  r@method <- "nlpca"
  r@missing <- sum(is.na(Matrix))

  R2cum <- rep(0, nPcs)
  for(i in 1:nPcs) 
    R2cum[i] <- 1 - sum((Matrix - fitted(r, Matrix, nPcs=i))^2, na.rm=TRUE) / sum(Matrix^2, na.rm=TRUE)

  R2 <- vector(length=length(R2cum), mode="numeric")
  R2[1] <- R2cum[1]
  if (length(R2cum) > 1) {
    for (i in 2:length(R2)) {
      R2[i] <- R2cum[i] - R2cum[i-1]
    }
  }
  r@R2 <- R2
  r@R2cum <- R2cum

  r
}



getHierarchicIdx <- function(hierarchicNum) {
  res <- matrix(1, ncol=hierarchicNum, nrow=hierarchicNum)
  res[lower.tri(res)] <- 0
  cbind(res, c(0, rep(1, hierarchicNum - 1)))
}

