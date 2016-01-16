##' Implements a Bayesian PCA missing value estimator.  The script
##' is a port of the Matlab version provided by Shigeyuki OBA.  See
##' also \url{http://ishiilab.jp/member/oba/tools/BPCAFill.html}.
##' BPCA combines an EM approach for PCA with a Bayesian model.  In
##' standard PCA data far from the training set but close to the
##' principal subspace may have the same reconstruction error.  BPCA
##' defines a likelihood function such that the likelihood for data
##' far from the training set is much lower, even if they are close to
##' the principal subspace.
##' 
##' Scores and loadings obtained with Bayesian PCA slightly differ
##' from those obtained with conventional PCA.  This is because BPCA
##' was developed especially for missing value estimation.  The
##' algorithm does not force orthogonality between factor loadings, as
##' a result factor loadings are not necessarily orthogonal.  However,
##' the BPCA authors found that including an orthogonality criterion
##' made the predictions worse.
##'
##' The authors also state that the difference between real and
##' predicted Eigenvalues becomes larger when the number of
##' observation is smaller, because it reflects the lack of
##' information to accurately determine true factor loadings from the
##' limited and noisy data.  As a result, weights of factors to
##' predict missing values are not the same as with conventional PCA,
##' but the missing value estimation is improved.
##' 
##' BPCA works iteratively, the complexity is growing with
##' \eqn{O(n^3)}{O(n^3)} because several matrix inversions are
##' required.  The size of the matrices to invert depends on the
##' number of components used for re-estimation.
##'
##' Finding the optimal number of components for estimation is not a
##' trivial task; the best choice depends on the internal structure of
##' the data.  A method called \code{kEstimate} is provided to
##' estimate the optimal number of components via cross validation.
##' In general few components are sufficient for reasonable estimation
##' accuracy. See also the package documentation for further
##' discussion about on what data PCA-based missing value estimation
##' makes sense.
##' 
##' It is not recommended to use this function directely but rather to
##' use the pca() wrapper function.
##'
##' There is a difference with respect the interpretation of rows
##' (observations) and columns (variables) compared to matlab
##' implementation. For estimation of missing values for microarray
##' data, the suggestion in the original bpca is to intepret genes as
##' observations and the samples as variables. In pcaMethods however,
##' genes are interpreted as variables and samples as observations
##' which arguably also is the more natural interpretation. For bpca
##' behavior like in the matlab implementation, simply transpose your
##' input matrix.
##' 
##' Details about the probabilistic model underlying BPCA are found in
##' Oba et. al 2003. The algorithm uses an expectation maximation
##' approach together with a Bayesian model to approximate the
##' principal axes (eigenvectors of the covariance matrix in PCA).
##' The estimation is done iteratively, the algorithm terminates if
##' either the maximum number of iterations was reached or if the
##' estimated increase in precision falls below \eqn{1e^{-4}}{1e^-4}.
##' 
##' \bold{Complexity:} The relatively high complexity of the method is
##' a result of several matrix inversions required in each step.
##' Considering the case that the maximum number of iteration steps is
##' needed, the approximate complexity is given by the term
##' \deqn{maxSteps \cdot row_{miss} \cdot O(n^3)}{maxSteps * row_miss
##' * O(n^3)} Where \eqn{row_{miss}}{row_miss} is the number of rows
##' containing missing values and \eqn{O(n^3)}{O(n^3)} is the
##' complexity for inverting a matrix of size
##' \eqn{components}{components}. Components is the number of
##' components used for re-estimation.
##' @title Bayesian PCA missing value estimation
##' @param Matrix \code{matrix} -- Pre-processed matrix (centered,
##'   scaled) with variables in columns and observations in rows. The
##'   data may contain missing values, denoted as \code{NA}.
##' @param nPcs \code{numeric} -- Number of components used for
##'   re-estimation. Choosing few components may decrease the
##'   estimation precision.
##' @param maxSteps \code{numeric} -- Maximum number of estimation
##'   steps.
##' @param verbose \code{boolean} -- BPCA prints the number of steps
##'   and the increase in precision if set to TRUE. Default is
##'   interactive().
##' @param threshold convergence threshold
##' @param ... Reserved for future use. Currently no further
##'   parameters are used
##' @return Standard PCA result object used by all PCA-based methods
##'   of this package. Contains scores, loadings, data mean and
##'   more. See \code{\link{pcaRes}} for details.
##' @references Shigeyuki Oba, Masa-aki Sato, Ichiro Takemasa, Morito
##'   Monden, Ken-ichi Matsubara and Shin Ishii.  A Bayesian missing
##'   value estimation method for gene expression profile
##'   data. \emph{Bioinformatics, 19(16):2088-2096, Nov 2003}.
##' @seealso \code{\link{ppca}}, \code{\link{svdImpute}},
##'   \code{\link{prcomp}}, \code{\link{nipalsPca}},
##'   \code{\link{pca}},
##'   \code{\link{pcaRes}}. \code{\link{kEstimate}}.
##' @note Requires \code{MASS}.
##' @examples
##' ## Load a sample metabolite dataset with 5\% missig values (metaboliteData)e
##' data(metaboliteData)
##' ## Perform Bayesian PCA with 2 components
##' pc <- pca(t(metaboliteData), method="bpca", nPcs=2)
##' ## Get the estimated principal axes (loadings)
##' loadings <- loadings(pc)
##' ## Get the estimated scores
##' scores <- scores(pc)
##' ## Get the estimated complete observations
##' cObs <- completeObs(pc)
##' ## Now make a scores and loadings plot
##' slplot(pc)
##' \dontshow{stopifnot(sum((fitted(pc) - t(metaboliteData))^2, na.rm=TRUE) < 200)}
##' @keywords multivariate
##' @export
##' @author Wolfram Stacklies
bpca <- function(Matrix, nPcs=2, maxSteps=100, 
                 verbose=interactive(), threshold=1e-4, ... ) {

  ## R implementation of a Bayesion PCA missing value estimator.
  ## After the Matlab script of Shigeyuki OBA (2002  May. 5th)
  ## See also: http://hawaii.aist-nara.ac.jp/%7Eshige-o/tools/
  ## Great thanks to them!
  M <- BPCA_initmodel(Matrix, nPcs)
  tauold <- 1000

  for( step in 1:maxSteps ) {
    M <- BPCA_dostep(M, Matrix)
    if( step %% 10 == 0 ) {
      tau <- M$tau
      dtau <- abs(log10(tau) - log10(tauold))
      if ( verbose ) {
        cat("Step Number           : ", step, '\n')
        cat("Increase in precision : ", dtau, '\n')
        cat("----------", '\n')
      }
      if (dtau < threshold) {
        break
      }
      tauold <- tau
    }
  }
  
  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for (i in 1:nPcs) {
    difference <-
      Matrix - (M$scores[,1:i, drop=FALSE] %*% t(M$PA[,1:i, drop=FALSE]) )
    R2cum[i] <- 1 - (sum(difference^2, na.rm=TRUE) / TSS)
  }

  result <- new("pcaRes")
  result@scores <- M$scores 
  result@loadings <- M$PA
  result@R2cum <- R2cum
  result@method <- "bpca"
  return(result)
}
