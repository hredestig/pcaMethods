\encoding{latin1}
\name{ppca}
\alias{ppca}
\title{Probabilistic PCA Missing Value Estimator}
\description{
  Implementation of probabilistic PCA (PPCA). PPCA allows to perform
  PCA on incomplete data and may be used for missing value
  estimation.
  This script was implemented after the Matlab version provided by
  Jakob Verbeek ( see \url{http://lear.inrialpes.fr/~verbeek/}) and
  the draft \emph{``EM Algorithms for PCA and Sensible PCA''} written
  by Sam Roweis. Thanks a lot! \cr

  Probabilistic PCA combines an EM approach for PCA with
  a probabilistic model. The EM approach is based on the assumption that
  the latent variables as well as the noise are normal distributed.

  In standard PCA data which is far from the training set but close to the
  principal subspace may have the same reconstruction error.
  PPCA defines a likelihood function such that the likelihood for data
  far from the training set is much lower, even if they are close to the
  principal subspace.
  This allows to improve the estimation accuracy.\cr

  A method called \code{kEstimate} is provided to estimate the optimal
  number of components via cross validation.
  In general few components are sufficient for reasonable estimation
  accuracy. See also the package documentation for further discussion
  on what kind of data PCA-based missing value estimation is advisable.\cr

  Requires \code{MASS}
}
\details{
  \bold{Complexity:} Runtime is linear in the number of data,
  number of data dimensions and number of principal components.\cr

  \bold{Convergence:} The algorithms seems not to converge to proper results
  in rare cases when unfavourable initial random numbers were chosen. 
  To avoid this you can set the seed (parameter seed) of the
  random number generator. \cr
  If used for missing value estimation, results may be checked by simply
  running the algorithm several times with changing seed, if the estimated
  values show little variance the algorithm converged well.
}
\usage{
  ppca(Matrix, nPcs, center = TRUE, completeObs = TRUE, seed = NA, ...)
}
\arguments{
  \item{Matrix}{\code{matrix} -- Data containing the variables in
    columns and observations in rows. The data may contain missing values,
    denoted as \code{NA}.}
  \item{nPcs}{\code{numeric} -- Number of components to estimate.
    The preciseness of the missing value estimation depends on the
    number of components, which should resemble the internal structure
    of the data.}
  \item{center}{\code{boolean} Mean center the data if TRUE}
  \item{completeObs}{\code{boolean} Return the complete observations if TRUE. This
    is the original data with NA values filled with the estimated values.}
  \item{seed}{\code{numeric} Set the seed for the random number generator. PPCA creates
    fills the initial loading matrix with random numbers chosen from a normal
    distribution. Thus results may vary slightly. Set the seed for
    exact reproduction of your results.}
  \item{...}{Reserved for future use. Currently no further parameters
    are used.}
}
\value{
  \item{pcaRes}{Standart PCA result object used by all
    PCA-based methods of this package. Contains scores, loadings, data mean and
    more. See \code{\link{pcaRes}} for details.}
}
\seealso{
  \code{\link{bpca}, \link{svdImpute}, \link{prcomp}, \link{nipalsPca}, \link{pca}, \link{pcaRes}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData)
data(metaboliteData)

# Now remove 10\% of the data
rows <- nrow(metaboliteData)
cols <- ncol(metaboliteData)
cond <- matrix(runif(rows * cols),rows,cols) < 0.1
metaboliteData[cond] <- NA

## Perform probabilistic PCA using the 3 largest components
result <- pca(metaboliteData, method="ppca", nPcs=3, center=TRUE)

## Get the estimated principal axes (loadings)
loadings <- result@loadings

## Get the estimated scores
scores <- result@scores

## Get the estimated complete observations
cObs <- result@completeObs

## Now plot the scores
plotPcs(result, scoresLoadings=c(TRUE,FALSE))

}
\keyword{multivariate}

\author{Wolfram Stacklies \cr
        Max Planck Institut fuer Molekulare Pflanzenphysiologie, Potsdam, Germany \cr
        \email{wolfram.stacklies@gmail.com} \cr
}

