##' This implements the SVDimpute algorithm as proposed by Troyanskaya
##' et al, 2001.  The idea behind the algorithm is to estimate the
##' missing values as a linear combination of the \code{k} most
##' significant eigengenes.
##'
##' Missing values are denoted as \code{NA}. It is not recommended
##' to use this function directely but rather to use the pca() wrapper
##' function.
##' 
##' As SVD can only be performed on complete matrices, all missing
##' values are initially replaced by 0 (what is in fact the mean on
##' centred data).  The algorithm works iteratively until the change
##' in the estimated solution falls below a certain threshold.  Each
##' step the eigengenes of the current estimate are calculated and
##' used to determine a new estimate. Eigengenes denote the loadings
##' if pca is performed considering variable (for Microarray data
##' genes) as observations.
##'
##' An optimal linear combination is found by regressing the
##' incomplete variable against the \code{k} most significant
##' eigengenes. If the value at position \code{j} is missing, the
##' \eqn{j^th}{j^th} value of the eigengenes is not used when
##' determining the regression coefficients.
##' @title SVDimpute algorithm
##' @param Matrix \code{matrix} -- Pre-processed (centered, scaled)
##' data with variables in columns and observations in rows. The data
##' may contain missing values, denoted as \code{NA}.
##' @param nPcs \code{numeric} -- Number of components to
##' estimate. The preciseness of the missing value estimation depends
##' on the number of components, which should resemble the internal
##' structure of the data.
##' @param threshold The iteration stops if the change in the matrix
##' falls below this threshold.
##' @param maxSteps Maximum number of iteration steps.
##' @param verbose Print some output if TRUE.
##' @param ... Reserved for parameters used in future version of the
##' algorithm
##' @note Each iteration, standard PCA (\code{prcomp}) needs to be
##' done for each incomplete variable to get the eigengenes. This is
##' usually fast for small data sets, but complexity may rise if the
##' data sets become very large.
##' @return Standard PCA result object used by all PCA-based methods
##' of this package. Contains scores, loadings, data mean and
##' more. See \code{\link{pcaRes}} for details.
##' @examples
##' ## Load a sample metabolite dataset with 5\% missing values
##' data(metaboliteData)
##' ## Perform svdImpute using the 3 largest components
##' result <- pca(metaboliteData, method="svdImpute", nPcs=3, center = TRUE)
##' ## Get the estimated complete observations
##' cObs <- completeObs(result)
##' ## Now plot the scores
##' plotPcs(result, type = "scores")
##' @keywords multivariate
##' @references Troyanskaya O. and Cantor M. and Sherlock G. and Brown
##' P. and Hastie T. and Tibshirani R. and Botstein D. and Altman
##' RB. - Missing value estimation methods for DNA
##' microarrays. \emph{Bioinformatics. 2001 Jun;17(6):520-5.}
##' @author Wolfram Stacklies
##' @export
svdImpute <- function(Matrix, nPcs=2, threshold=0.01, maxSteps=100,
                      verbose=interactive(), ...) {
  missing <- is.na(Matrix)
  temp <- apply(missing, 2, sum)
  missIx <- which(temp != 0)

  ## Initially set estimates to 0
  Matrix[missing] <- 0

  ## Now do the regression
  count <- 0
  error <- Inf

  while ( (error > threshold) && (count < maxSteps) ) {
    res <- prcomp(t(Matrix), center = FALSE, scale = FALSE, retx = TRUE)
    loadings <- res$rotation[,1:nPcs, drop = FALSE]
    sDev <- res$sdev

    ## Estimate missing values as a linear combination of the eigenvectors
    ## The optimal solution is found by regression against the k eigengenes
    for (index in missIx) {
      target <- Matrix[!missing[,index],index, drop = FALSE]
      Apart <- loadings[!missing[,index], , drop = FALSE]
      Bpart <- loadings[missing[,index], , drop = FALSE]
      X <- MASS::ginv(Apart) %*% target
      estimate <- Bpart %*% X
      
      Matrix[missing[,index], index] <- estimate
    }
    
    count <- count + 1
    if (count > 5) {
      error <- sqrt(sum( (MatrixOld - Matrix)^2 ) / sum(MatrixOld^2))
      if (verbose) { cat("change in estimate: ", error, "\n") }
    }
    MatrixOld <- Matrix
  }

  tmp <- prcomp(Matrix, center = FALSE, scale = FALSE, retx = TRUE)
  loadings <- cbind(tmp$rotation[,1:nPcs])
  scores <- cbind(tmp$x[,1:nPcs])

  ## Calculate R2cum
  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for (i in 1:nPcs) {
    difference <-
      Matrix - (scores[,1:i, drop=FALSE] %*% t(loadings[,1:i, drop=FALSE]))
    R2cum[i] <- 1 - (sum(difference^2) / TSS)
  }

  result <- new("pcaRes")
  result@scores <- scores
  result@loadings <- loadings
  result@R2cum <- R2cum
  result@method <- "svdImpute"

  return(result)
}
