##' This is a PCA implementation robust to outliers in a data set. It
##' can also handle missing values, it is however NOT intended to be
##' used for missing value estimation.  As it is based on robustSVD we
##' will get an accurate estimation for the loadings also for
##' incomplete data or for data with outliers.  The returned scores
##' are, however, affected by the outliers as they are calculated
##' inputData X loadings. This also implies that you should look at
##' the returned R2/R2cum values with caution.  If the data show
##' missing values, scores are caluclated by just setting all NA -
##' values to zero. This is not expected to produce accurate results.
##' Please have also a look at the manual page for \code{robustSvd}.
##' Thus this method should mainly be seen as an attempt to integrate
##' \code{robustSvd()} into the framework of this package.  Use one of
##' the other methods coming with this package (like PPCA or BPCA) if
##' you want to do missing value estimation.  It is not recommended to
##' use this function directely but rather to use the pca() wrapper
##' function.
##' 
##' The method is very similar to the standard \code{prcomp()}
##' function.  The main difference is that \code{robustSvd()} is used
##' instead of the conventional \code{svd()} method.
##' @title PCA implementation based on robustSvd
##' @param Matrix \code{matrix} -- Data containing the variables in
##' columns and observations in rows. The data may contain missing
##' values, denoted as \code{NA}.
##' @param nPcs \code{numeric} -- Number of components to
##' estimate. The preciseness of the missing value estimation depends
##' on the number of components, which should resemble the internal
##' structure of the data.
##' @param verbose \code{boolean} Print some output to the command
##' line if TRUE
##' @param ...  Reserved for future use. Currently no further
##' parameters are used
##' @return Standard PCA result object used by all PCA-based methods
##' of this package. Contains scores, loadings, data mean and
##' more. See \code{\link{pcaRes}} for details.  are used.
##' @seealso \code{\link{robustSvd}, \link{svd}, \link{prcomp},
##' \link{pcaRes}}.
##' @examples
##' ## Load a complete sample metabolite data set and mean center the data
##' data(metaboliteDataComplete)
##' mdc <- scale(metaboliteDataComplete, center=TRUE, scale=FALSE)
##' ## Now create 5\% of outliers.
##' cond   <- runif(length(mdc)) < 0.05;
##' mdcOut <- mdc
##' mdcOut[cond] <- 10
##' ## Now we do a conventional PCA and robustPca on the original and the data
##' ## with outliers.
##' ## We use center=FALSE here because the large artificial outliers would
##' ## affect the means and not allow to objectively compare the results.
##' resSvd    <- pca(mdc, method="svd", nPcs=10, center=FALSE)
##' resSvdOut <- pca(mdcOut, method="svd", nPcs=10, center=FALSE)
##' resRobPca <- pca(mdcOut, method="robustPca", nPcs=10, center=FALSE)
##' ## Now we plot the results for the original data against those with outliers
##' ## We can see that robustPca is hardly effected by the outliers.
##' plot(loadings(resSvd)[,1], loadings(resSvdOut)[,1])
##' plot(loadings(resSvd)[,1], loadings(resRobPca)[,1])
##' @keywords multivariate
##' @export
##' @author Wolfram Stacklies
robustPca <- function(Matrix, nPcs=2, verbose=interactive(), ... ) {

  nas <- is.na(Matrix)

  if (sum(nas) != 0)
    warning("Data is incomplete, it is not recommended to use robustPca for missing value estimation")
  svdSol <- robustSvd(Matrix)
  
  ## Sort the eigenvalues and eigenvectors
  loadings <- svdSol$v[, 1:nPcs, drop=FALSE]
  sDev     <- svdSol$d[1:nPcs] / sqrt(max(1, nrow(Matrix) - 1))

  ## We estimate the scores by just setting all NA values to 0 This is
  ## a bad approximation, I know... Use ppca / bpca or other missing
  ## value estimation methods included in this package
  compMat <- Matrix
  compMat[is.na(compMat)] <- 0
  scores   <- compMat %*% loadings

  ## Calculate R2cum (on the complete observations only)
  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for (i in 1:nPcs) {
    difference <- Matrix -
      (scores[,1:i, drop=FALSE] %*% t(loadings[,1:i, drop=FALSE]))
    R2cum[i] <- 1 - (sum(difference^2) / TSS)
  }

  result <- new("pcaRes")
  result@loadings <- loadings
  result@scores <- scores
  result@R2cum <- R2cum
  result@method <- "robustPca"
  return(result)
}

##' A robust approximation to the singular value decomposition of a
##' rectangular matrix is computed using an alternating L1 norm
##' (instead of the more usual least squares L2 norm). As the SVD is
##' a least-squares procedure, it is highly susceptible to outliers
##' and in the extreme case, an individual cell (if sufficiently
##' outlying) can draw even the leading principal component toward
##' itself.
##'
##' See Hawkins et al (2001) for details on the robust SVD algorithm.
##' Briefly, the idea is to sequentially estimate the left and right
##' eigenvectors using an L1 (absolute value) norm minimization.
##'
##' Note that the robust SVD is able to accomodate missing values in
##' the matrix \code{x}, unlike the usual \code{svd} function.
##' 
##' Also note that the eigenvectors returned by the robust SVD
##' algorithm are NOT (in general) orthogonal and the eigenvalues need
##' not be descending in order.
##' @title Alternating L1 Singular Value Decomposition
##' @param x A matrix whose SVD decomposition is to be
##' computed. Missing values are allowed.
##' @return The robust SVD of the matrix is x=u d v'. \item{d}{A
##' vector containing the singular values of \code{x}.} \item{u}{A
##' matrix whose columns are the left singular vectors of \code{x}.}
##' \item{v}{A matrix whose columns are the right singular vectors of
##' \code{x}.}
##' @note Two differences from the usual SVD may be noted. One relates
##' to orthogonality. In the conventional SVD, all the eigenvectors
##' are orthogonal even if not explicitly imposed.  Those returned by
##' the AL1 algorithm (used here) are (in general) not orthogonal.
##' Another difference is that, in the L2 analysis of the conventional
##' SVD, the successive eigen triples (eigenvalue, left eigenvector,
##' right eigenvector) are found in descending order of
##' eigenvalue. This is not necessarily the case with the AL1
##' algorithm.  Hawkins et al (2001) note that a larger eigen value
##' may follow a smaller one.
##' @references Hawkins, Douglas M, Li Liu, and S Stanley Young (2001)
##' Robust Singular Value Decomposition, National Institute of
##' Statistical Sciences, Technical Report Number
##' 122. \url{http://www.niss.org/technicalreports/tr122.pdf}
##' @author Kevin Wright, modifications by Wolfram Stacklies
##' @seealso \code{\link{svd}}, \code{\link[ade4:nipals]{nipals}} for
##' an alternating L2 norm method that also accommodates missing data.
##' @examples
##' ## Load a complete sample metabolite data set and mean center the data
##' data(metaboliteDataComplete)
##' mdc <- prep(metaboliteDataComplete, center=TRUE, scale="none")
##' ## Now create 5% of outliers.
##' cond   <- runif(length(mdc)) < 0.05;
##' mdcOut <- mdc
##' mdcOut[cond] <- 10
##' ## Now we do a conventional SVD and a robustSvd on both, the original and the 
##' ## data with outliers.
##' resSvd <- svd(mdc)
##' resSvdOut <- svd(mdcOut)
##' resRobSvd <- robustSvd(mdc)
##' resRobSvdOut <- robustSvd(mdcOut)
##' ## Now we plot the results for the original data against those with outliers
##' ## We can see that robustSvd is hardly affected by the outliers.
##' plot(resSvd$v[,1], resSvdOut$v[,1])
##' plot(resRobSvd$v[,1], resRobSvdOut$v[,1])
##' @keywords algebra
##' @export
robustSvd <- function(x) {
  ## We need the weightedMedian function provided by the aroma.light
  ## package. However we do not want to make the whole package dependant
  ## on aroma.light
  if (!requireNamespace("matrixStats", quietly=TRUE))
    stop("package matrixStats required but not available")

  L1RegCoef <- function(x,a){
    keep <- (a!=0) & (!is.na(x))
    a <- a[keep]
    return(matrixStats::weightedMedian(x[keep]/a, abs(a),
                                       na.rm=TRUE, interpolate=FALSE))
  }

  L1Eigen <- function(x,a,b){
    x <- as.vector(x) # Convert from matrix to vector
    ab <- as.vector(outer(a,b))
    keep <- (ab!=0) & (!is.na(x))
    ab <- ab[keep]
    return(matrixStats::weightedMedian(x[keep]/ab, abs(ab),
                                       na.rm=TRUE, interpolate=FALSE))
  }

  ## Initialize outputs
  svdu <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  svdv <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
  svdd <- rep(NA,ncol(x))

  for(k in 1:ncol(x)) {
    ak <- apply(abs(x),1,median,na.rm=TRUE)
    converged <- FALSE
    
    while(!converged) {
      akprev <- ak
      c <- apply(x,2,L1RegCoef,ak)
      bk <- c/sqrt(sum(c^2))
      d <- apply(x,1,L1RegCoef,bk)
      ak <- d/sqrt(sum(d^2))
      if(sum((ak-akprev)^2)< 1e-10) converged <- TRUE
    }
    eigenk <- L1Eigen(x,ak,bk)
    ## Deflate the x matrix
    x <- x - eigenk * ak %*% t(bk)
    ## Store eigen triple for output
    svdu[,k] <- ak
    svdv[,k] <- bk 
    svdd[k] <- eigenk
  }

  ## Create the result object
  ret <- list()

  ret$d <- svdd
  ret$u <- svdu
  ret$v <- svdv
  return(ret)
}

