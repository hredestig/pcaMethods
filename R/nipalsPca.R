##' PCA by non-linear iterative partial least squares
##'
##' Can be used for computing PCA on a numeric matrix using either the
##' NIPALS algorithm which is an iterative approach for estimating the
##' principal components extracting them one at a time. NIPALS can
##' handle a small amount of missing values. It is not recommended to
##' use this function directely but rather to use the pca() wrapper
##' function. There is a C++ implementation given as \code{nipalsPca} which is faster.
##' @title NIPALS PCA implemented in R
##' @param Matrix Numerical matrix samples in rows and variables as columns.
##' @param nPcs Number of components that should be extracted.
##' @param center Mean center the data column wise if set TRUE
##' @param completeObs Return the estimated complete observations. This is
##' the input Matrix with NA values replaced by the estimated values.
##' @param varLimit Optionally the ratio of variance that should be
##' explained. \code{nPcs} is ignored if varLimit < 1
##' @param maxSteps Defines how many iterations can be done before
##' algorithm should abort (happens almost exclusively when there were
##' some wrong in the input data).
##' @param threshold The limit condition for judging if the algorithm has
##' converged or not, specifically if a new iteration is done if
##' \eqn{(T_{old} - T)^T(T_{old} - T) > \code{limit}}.
##' @param verbose Show simple progress information.
##' @param ... Only used for passing through arguments.
##' @return A \code{pcaRes} object.
##' @references
##' Wold, H. (1966) Estimation of principal components and related models by
##' iterative least squares. In Multivariate Analysis (Ed.,
##' P.R. Krishnaiah), Academic Press, NY, 391-420.
##' @author Henning Redestig
##' @seealso \code{prcomp}, \code{princomp}, \code{pca}
##' @examples
##' data(iris)
##' mat <- as.matrix(iris[,1:4])
##' system.time(pcIr <- RnipalsPca(mat, method="rnipals", nPcs=3))
##' system.time(pcIr <- nipalsPca(mat, nPcs=3))
##' @keywords multivariate
RnipalsPca <- function(Matrix, nPcs=2, center=TRUE, completeObs=TRUE,
                      varLimit=1,
                      maxSteps=5000, 
                      threshold=1e-6, verbose=interactive(), ...) {

  if (center) {
    object <- scale(Matrix, center = TRUE, scale = FALSE)
    means <- attr(object, "scaled:center")
  } else
  object <- Matrix

  missing <- is.na(Matrix)
  nObs <- nrow(object)
  nVar <- ncol(object)

  ##Find a good? starting column -- better way?
  startingColumn <- 1

  ## sum(c(NA, NA), na.rm=TRUE) is 0, but we want NA
  sum.na <- function(x){ ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))}

  TotalSS <- sum(object*object, na.rm=TRUE)

  ph <- rep(0, nVar)
  R2cum <- NULL
  scores <- NULL
  loadings <- NULL
  anotherPc <- TRUE
  l <- 1
  
  while(anotherPc) {
    count <- 0                 #number of iterations done
    th <- object[,startingColumn]   #first column is starting vector for th
    continue <- TRUE
    if(verbose) cat(paste("Calculating PC", l, ": ", sep=""))
    
    while(continue) {
      count <- count+1
      ph <- rep(0, nVar)

      ##Calculate loadings through LS regression
      ##Note: object*th is column-wise multiplication
      tsize <- sum(th * th, na.rm=TRUE)
      ph <- apply(object * (th / tsize), 2, sum.na)
      ##normalize ph based on the available values.
      psize <- sum(ph*ph, na.rm=TRUE)
      ph <- ph / sqrt(psize)
      
      ##Calculate scores through LS regression
      ##Trick: To get row-wise multiplication, use t(object)*ph, then
      ##be sure to use apply(,2,) and NOT apply(,1,)!
      th.old <- th
      th <- apply(t(object) * ph, 2, sum.na)
      
      ##Round up by calculating if convergence condition is met and
      ##checking if it seems to be an neverending loop.
      if (count > maxSteps) {
        stop("Too many iterations, quitting")
      }
      if (t(na.omit(th.old - th)) %*% (na.omit(th.old - th)) <= threshold) {
        continue = FALSE
      }
      if (verbose)cat("*")
    }
    if (verbose) cat(" Done\n")
    object <- object - (th %*% t(ph))
    scores <- cbind(scores, th)
    loadings <- cbind(loadings, ph)
    
    ##cumulative proportion of variance
    R2cum <- cbind(R2cum, 1 - (sum(object*object,na.rm=TRUE) / TotalSS))
    l <- l + 1
    if((!abs(varLimit - 1) < 1e-4 & R2cum[1,l - 1] >= varLimit)| l > nPcs) {
      anotherPc <- FALSE
      nPcs <- l - 1
    }
  }
  R2 <- R2cum[1]
  if(length(R2cum) > 1) R2 <- c(R2, diff(R2cum[1,]))

  if (completeObs) {
    Ye <- scores %*% t(loadings)
    if (center) {
      Ye <- Ye + rep(means, each=nrow(Ye)) # Addition is column-wise
    }
    cObs <- Matrix
    cObs[missing] <- Ye[missing]
  }

  rownames(scores) <- rownames(object)
  colnames(scores) <- paste("PC", 1:nPcs, sep="")
  rownames(loadings) <- colnames(object)
  colnames(loadings) <- paste("PC", 1:nPcs, sep="")
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- cObs
  r@scores <- scores
  r@loadings <- loadings
  r@R2cum <- c(R2cum)
  r@sDev <- apply(scores, 2, sd)
  r@R2 <- R2
  r@nObs <- nObs
  r@nVar <- nVar
  r@varLimit <- varLimit
  r@nPcs <- nPcs
  r@centered <- center
  r@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
  r@method <- "rnipals"
  r@missing <- sum(is.na(Matrix))
  return(r)
}


##' PCA by non-linear iterative partial least squares
##'
##' Can be used for computing PCA on a numeric matrix using either the
##' NIPALS algorithm which is an iterative approach for estimating the
##' principal components extracting them one at a time. NIPALS can
##' handle a small amount of missing values. It is not recommended to
##' use this function directely but rather to use the pca() wrapper
##' function.
##' @title NIPALS PCA
##' @param Matrix Numerical matrix samples in rows and variables as columns.
##' @param nPcs Number of components that should be extracted.
##' @param center Mean center the data column wise if set TRUE
##' @param completeObs Return the estimated complete observations. This is
##' the input Matrix with NA values replaced by the estimated values.
##' @param varLimit Optionally the ratio of variance that should be
##' explained. \code{nPcs} is ignored if varLimit < 1
##' @param maxSteps Defines how many iterations can be done before
##' algorithm should abort (happens almost exclusively when there were
##' some wrong in the input data).
##' @param threshold The limit condition for judging if the algorithm has
##' converged or not, specifically if a new iteration is done if
##' \eqn{(T_{old} - T)^T(T_{old} - T) > \code{limit}}.
##' @param ... Only used for passing through arguments.
##' @return A \code{pcaRes} object.
##' @references
##' Wold, H. (1966) Estimation of principal components and related models by
##' iterative least squares. In Multivariate Analysis (Ed.,
##' P.R. Krishnaiah), Academic Press, NY, 391-420.
##' @author Henning Redestig
##' @seealso \code{prcomp}, \code{princomp}, \code{pca}
##' @examples
##' data(iris)
##' mat <- as.matrix(iris[,1:4])
##' pcIr <- nipalsPca(mat, nPcs=2)
##' @keywords multivariate
nipalsPca <- function(Matrix, nPcs=2, center=TRUE, completeObs=TRUE,
                      varLimit=1,
                      maxSteps=5000, 
                      threshold=1e-6, ...) {

  if (center) {
    object <- scale(Matrix, center = TRUE, scale = FALSE)
    means <- attr(object, "scaled:center")
  } else
  object <- Matrix

  missing <- is.na(Matrix)
  nObs <- nrow(object)
  nVar <- ncol(object)
  
  nipRes <- .Call("Nipals", Matrix,
                  params=list(nPcs=nPcs,
                    varLimit=varLimit,
                    threshold=threshold,
                    maxSteps=maxSteps), PACKAGE="pcaMethods")

  R2 <- nipRes$R2
  scores <- nipRes$scores
  loadings <- nipRes$loadings

  R2cum <- cumsum(R2)

  if (completeObs) {
    Ye <- scores %*% t(loadings)
    if (center)
      Ye <- sweep(Ye, 2, means, "+")
    cObs <- Matrix
    cObs[missing] <- Ye[missing]
  }

  rownames(scores) <- rownames(object)
  colnames(scores) <- paste("PC", 1:nPcs, sep="")
  rownames(loadings) <- colnames(object)
  colnames(loadings) <- paste("PC", 1:nPcs, sep="")
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- cObs
  r@scores <- scores
  r@loadings <- loadings
  r@R2cum <- c(R2cum)
  r@sDev <- apply(scores, 2, sd)
  r@R2 <- R2
  r@nObs <- nObs
  r@nVar <- nVar
  r@varLimit <- varLimit
  r@nPcs <- nPcs
  r@centered <- center
  r@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
  r@method <- "nipals"
  r@missing <- sum(is.na(Matrix))
  return(r)
}


