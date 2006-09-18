\encoding{latin1}
\name{bpca}
\alias{bpca}
\title{Bayesian PCA Missing Value Estimator}
\description{
	Implements a Bayesian PCA missing value estimator.
	The script is a port of the Matlab version provided by
	Shigeyuki OBA.
	See also \url{http://hawaii.aist-nara.ac.jp/\%7Eshige-o/tools/}.\cr
	BPCA combines an EM approach for PCA with
	a Bayesian model.
	In standard PCA data far from the training set but close to the
	principal subspace may have the same reconstruction error.
	BPCA defines a likelihood function such that the likelihood for data
	far from the training set is much lower, even if they are close to the
	principal subspace.
	This allows to improve the estimation accuracy.

	The method works iteratively. Its complexity is growing with
	\eqn{O(n^3)}{O(n^3)} because several matrix inversions are required.
	The size of the matrices to invert depends on the number of components
	used for re-estimation.\cr
	Finding the optimal number of components for estimation is not a
	trivial task; the best choice depends on the internal structure of the
	data.
	A method called \code{kEstimate} is provided to estimate the optimal
	number of components via cross validation.
	In general few components are sufficient for reasonable estimation
	accuracy. See also the package documentation for further discussion
	on what kind of data PCA-based missing value estimation is advisable.\cr
	Please not that results obtained with \code{Q2} cross validation may
	not reflect the method's actual accuracy because the signs of scores
	and loadings may be different for independent runs.\cr

	Requires \code{MASS}.
}
\details{
	Details about the probabilistic model underlying BPCA are found in
	Oba et. al 2003. The algorithm uses an expectation maximation approach
	together with a Bayesian model to approximate the principal axes
	(eigenvectors of the covariance matrix in PCA).
	The estimation is done iteratively, the algorithm terminates if either
	the maximum number of iterations was reached or if the estimated
	increase in precision falls below \eqn{1e^{-4}}{1e^-4}.

	\bold{Complexity:} The relatively high complexity of the method is a result
	of several matrix inversions required in each step.
	Considering the case that the maximum number of iteration steps is needed,
	the approximate complexity is given by the term
	\deqn{maxSteps \cdot row_{miss} \cdot O(n^3)}{maxSteps * row_miss * O(n^3)}
	Where \eqn{row_{miss}}{row_miss} is the number of rows containing missing values
	and \eqn{O(n^3)}{O(n^3)} is the complexity for inverting a matrix of size
	\eqn{components ×components}. Components is the number of components
	used for re-estimation.
}
\usage{
        bpca(Matrix, nPcs = 2, maxSteps = 100, verbose = interactive(), ...)
}
\arguments{
	\item{Matrix}{\code{matrix} -- Data containing the variables in
	columns and observations in rows. The data may contain missing values,
	denoted as \code{NA}.}
	\item{nPcs}{\code{numeric} -- Number of components used for re-estimation.
	Choosing few components may decrease the estimation
	precision.}
	\item{maxSteps}{\code{numeric} -- Maximum number of estimation steps.
	Default is 100. }
	\item{verbose}{\code{boolean} -- BPCA prints the number of steps and the
	increase in precision if set to TRUE. Default is interactive().}
	\item{...}{Reserved for future use. Currently no further parameters are
	used}
}
\value{
	\item{pcaRes}{Standard PCA result object used by all
	PCA-based methods of this package. Contains scores, loadings, data mean and
	more. See \code{\link{pcaRes}} for details.}
}
\references{
	Shigeyuki Oba, Masa-aki Sato, Ichiro Takemasa,
	Morito Monden, Ken-ichi Matsubara and Shin Ishii.
	A Bayesian missing value estimation method for gene expression profile
	data. \emph{Bioinformatics, 19(16):2088-2096, Nov 2003}.
}
\seealso{
  \code{\link{ppca}, \link{svdImpute}, \link{prcomp}, \link{nipalsPca}, \link{pca}, \link{pcaRes}. \link{kEstimate}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData)
data(metaboliteData)

# Now remove 10\% of the data
rows <- nrow(metaboliteData)
cols <- ncol(metaboliteData)
cond<-matrix(runif(rows * cols),rows,cols) < 0.1
metaboliteData[cond] <- NA

## Perform Bayesian PCA with 3 components
result <- pca(metaboliteData, method="bpca", nPcs=2, center=FALSE, scale="none")

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
