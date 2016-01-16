##' Implementation of probabilistic PCA (PPCA). PPCA allows to perform
##' PCA on incomplete data and may be used for missing value
##' estimation.  This script was implemented after the Matlab version
##' provided by Jakob Verbeek ( see
##' \url{http://lear.inrialpes.fr/~verbeek/}) and the draft \emph{``EM
##' Algorithms for PCA and Sensible PCA''} written by Sam Roweis.
##'
##' Probabilistic PCA combines an EM approach for PCA with a
##' probabilistic model. The EM approach is based on the assumption
##' that the latent variables as well as the noise are normal
##' distributed.
##'
##' In standard PCA data which is far from the training set but close
##' to the principal subspace may have the same reconstruction error.
##' PPCA defines a likelihood function such that the likelihood for
##' data far from the training set is much lower, even if they are
##' close to the principal subspace.  This allows to improve the
##' estimation accuracy.
##'
##' A method called \code{kEstimate} is provided to estimate the
##' optimal number of components via cross validation.  In general few
##' components are sufficient for reasonable estimation accuracy. See
##' also the package documentation for further discussion on what kind
##' of data PCA-based missing value estimation is advisable.
##'
##' \bold{Complexity:}\cr Runtime is linear in the number of data,
##' number of data dimensions and number of principal components.
##'
##' \bold{Convergence:}  The threshold indicating convergence was
##' changed from 1e-3 in 1.2.x to 1e-5 in the current version  leading
##' to more stable results.  For reproducability you can set the seed
##' (parameter seed) of the random number generator. If used for
##' missing value estimation, results may be checked by simply running
##' the algorithm several times with changing seed, if the estimated
##' values show little variance the algorithm converged well. 
##' @title Probabilistic PCA 
##' @param Matrix \code{matrix} -- Data containing the variables in
##' columns and observations in rows. The data may contain missing
##' values, denoted as \code{NA}.
##' @param nPcs \code{numeric} -- Number of components to
##' estimate. The preciseness of the missing value estimation depends
##' on the number of components, which should resemble the internal
##' structure of the data.
##' @param seed \code{numeric} Set the seed for the random number
##' generator. PPCA creates fills the initial loading matrix with
##' random numbers chosen from a normal distribution. Thus results may
##' vary slightly. Set the seed for exact reproduction of your
##' results.
##' @param threshold Convergence threshold.
##' @param maxIterations the maximum number of allowed iterations
##' @param ... Reserved for future use. Currently no further
##' parameters are used.
##' @note Requires \code{MASS}. It is not recommended to use this
##' function directely but rather to use the pca() wrapper function.
##' @return Standard PCA result object used by all PCA-based methods
##' of this package. Contains scores, loadings, data mean and
##' more. See \code{\link{pcaRes}} for details.
##' @seealso \code{\link{bpca}, \link{svdImpute}, \link{prcomp},
##' \link{nipalsPca}, \link{pca}, \link{pcaRes}}.
##' @examples
##' ## Load a sample metabolite dataset with 5\% missing values (metaboliteData)
##' data(metaboliteData)
##' ## Perform probabilistic PCA using the 3 largest components
##' result <- pca(t(metaboliteData), method="ppca", nPcs=3, seed=123)
##' ## Get the estimated complete observations
##' cObs <- completeObs(result)
##' ## Plot the scores
##' plotPcs(result, type = "scores")
##' \dontshow{
##'   stopifnot(sum((fitted(result) - t(metaboliteData))^2, na.rm=TRUE) < 200)
##' }
##' @keywords multivariate
##' @author Wolfram Stacklies
##' @export
ppca <- function(Matrix, nPcs=2, seed=NA, threshold=1e-5, maxIterations=1000, ...) {
  ## Set the seed to the user defined value. This affects the generation
  ## of random values for the initial setup of the loading matrix
  if (!is.na(seed)) 
    set.seed(seed)

  N <- nrow(Matrix)
  D <- ncol(Matrix)

  Obs <- !is.na(Matrix)
  hidden <- which(is.na(Matrix))
  missing <- length(hidden)

  if(missing) { Matrix[hidden] <- 0 } 

  ## ------- Initialization
  r <- sample(N)
  C <- t(Matrix[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as Matrix
  C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  CtC <- t(C) %*% C
  ## inv(C'C) C' X is the solution to the EM problem
  X  <- Matrix %*% C %*% solve(CtC)
  recon    <- X %*% t(C)
  recon[hidden] <- 0
  ss <- sum(sum((recon - Matrix)^2)) / (N * D - missing)

  count <- 1
  old <- Inf
  
  ## ------ EM iterations
  while (count > 0) {
    ## E-step, (co)variances
    Sx <- solve(diag(nPcs) + CtC/ss) 
    ss_old <- ss
    if(missing) {
      proj <- X %*% t(C)
      Matrix[hidden] <- proj[hidden]
    }
    
    ## E step: expected values
    X <- Matrix %*% C %*% Sx / ss
    
    ## M-step
    SumXtX <- t(X) %*% X

    ## Replace the right matrix division from matlab
    C <- (t(Matrix) %*% X) %*% solve( (SumXtX + N * Sx) )
    
    CtC <- t(C) %*% C
    ss <- ( sum(sum( (C %*% t(X) - t(Matrix))^2 )) + N * sum(sum(CtC %*% Sx)) +
           missing * ss_old ) / (N * D)

    objective <- N * (D * log(ss) + sum(diag(Sx)) - log(det(Sx)) ) +
      sum(diag(SumXtX)) - missing * log(ss_old)

    rel_ch <- abs( 1 - objective / old )
    old <- objective

    count <- count + 1
    if( rel_ch < threshold & count > 5 ) {
      count <- 0
    }
    else if (count > maxIterations) {
      count <- 0
      warning("stopped after max iterations, but rel_ch was > threshold")
    }
  } ## End EM iteration
  C <- orth(C)
  evs <- eigen( cov(Matrix %*% C) )
  vals <- evs[[1]]
  vecs <- evs[[2]]
  
  C <- C %*% vecs
  X <- Matrix %*% C

  ## Paramters in original Matlab implementation were:
  ## C (D by d)    - C has the approximate loadings (eigenvectors of
  ## the covariance matrix)
  ##          as columns.
  ## X        - The approximate scores 
  ## Matrix (N by D)    - Expected complete observations.
  ## M (D by 1)    - Column wise data mean
  ## ss (scalar)    - isotropic variance outside subspace

  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for (i in 1:ncol(C)) {
    difference <- Matrix - (X[,1:i, drop=FALSE] %*% t(C[,1:i, drop=FALSE]))
    R2cum[i] <- 1 - (sum(difference^2, na.rm=TRUE) / TSS)
  }

  res <- new("pcaRes")
  res@scores <- X
  res@loadings <- C
  res@R2cum <- R2cum
  res@method <- "ppca"
  return(res)
}
