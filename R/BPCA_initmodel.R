##' Model initialization for Bayesian PCA. This function is NOT
##' inteded to be run separately!
##'
##' The function calculates the initial Eigenvectors by use of SVD
##' from the complete rows.  The data structure M is created and
##' initial values are  assigned.
##' @title Initialize BPCA model
##' @param y numeric matrix containing missing values. Missing values
##' are denoted as 'NA'
##' @param components Number of components used for estimation
##' @return List containing
##' \item{rows}{Row number of input matrix}
##' \item{cols}{Column number of input matrix}
##' \item{comps}{Number of components to use}
##' \item{yest}{(working variable) current estimate of complete data}
##' \item{row_miss}{(Array) Indizes of rows containing missing values}
##' \item{row_nomiss}{(Array) Indices of complete rows (such with no
##' missing values)}
##' \item{nans}{Matrix of same size as input data. TRUE if \code{input == NA},
##' false otherwise}
##' \item{mean}{Column wise data mean}
##' \item{PA}{ (d x k) Estimated principal axes (eigenvectors,
##' loadings) The matrix ROWS are the vectors}
##' \item{tau}{Estimated precision of the residual error}
##' \item{scores}{ Estimated scores}
##' Further elements are: galpha0, balpha0, alpha, gmu0, btau0, gtau0,
##' SigW. These are working variables or constants.
##' @author Wolfram Stacklies
BPCA_initmodel <- function(y, components) {
  ## Initialization, write static parameters to the central
  M <- NULL 
  M$rows <- nrow(y)
  M$cols <- ncol(y) 
  M$comps <- components ## Column number
  M$yest <- y ## Original data, NAs are set to 0 later on

  ## Find rows with missing values, etc...
  M$nans <- is.na(y)
  temp <- apply(M$nans, 1, sum)
  M$row_nomiss <- which(temp == 0)
  M$row_miss <- which(temp != 0)
  M$yest[M$nans] <- 0
  M$scores <- NULL

  ## Get the SVD of the complete rows
  covy <- cov(M$yest)
  values <- svd(covy, components, components)
  U <- values[[2]]
  S <- diag( values[[1]][1:components], nrow = components, ncol = components)
  V <- values[[3]]

  ## M$mean: column wise mean of the original data
  M$mean <- matrix(0, 1, M$cols)
  for(j in 1:M$cols) {
    idx <- which(!is.na(y[,j]))
    M$mean[j] <- mean(y[idx,j])
  }

  M$PA <- U %*% sqrt(S)
  M$tau <- 1 / ( sum(diag(covy)) - sum(diag(S)) )
  
  ## Constants etc
  taumax <- 1e10
  taumin <- 1e-10
  M$tau <- max( min(M$tau, taumax), taumin )

  M$galpha0 <- 1e-10
  M$balpha0 <- 1
  M$alpha <- (2 * M$galpha0 + M$cols) / (M$tau * diag(t(M$PA) %*% M$PA) + 2 * M$galpha0 / M$balpha0)

  M$gmu0 <- 0.001

  M$btau0 <- 1
  M$gtau0 <- 1e-10
  M$SigW <- diag(components)
  return(M)
}
