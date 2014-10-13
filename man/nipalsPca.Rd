\name{nipalsPca}
\alias{nipalsPca}
\title{NIPALS PCA}
\usage{
  nipalsPca(Matrix, nPcs = 2, varLimit = 1,
    maxSteps = 5000, threshold = 1e-06, ...)
}
\arguments{
  \item{Matrix}{Pre-processed (centered, scaled) numerical
  matrix samples in rows and variables as columns.}

  \item{nPcs}{Number of components that should be
  extracted.}

  \item{varLimit}{Optionally the ratio of variance that
  should be explained. \code{nPcs} is ignored if varLimit <
  1}

  \item{maxSteps}{Defines how many iterations can be done
  before algorithm should abort (happens almost exclusively
  when there were some wrong in the input data).}

  \item{threshold}{The limit condition for judging if the
  algorithm has converged or not, specifically if a new
  iteration is done if \eqn{(T_{old} - T)^T(T_{old} - T) >
  \code{limit}}.}

  \item{...}{Only used for passing through arguments.}
}
\value{
  A \code{pcaRes} object.
}
\description{
  PCA by non-linear iterative partial least squares
}
\details{
  Can be used for computing PCA on a numeric matrix using
  either the NIPALS algorithm which is an iterative
  approach for estimating the principal components
  extracting them one at a time. NIPALS can handle a small
  amount of missing values. It is not recommended to use
  this function directely but rather to use the pca()
  wrapper function.
}
\examples{
data(metaboliteData)
mat <- prep(t(metaboliteData))
pc <- nipalsPca(mat, nPcs=2)
## better use pca()
pc <- pca(t(metaboliteData), method="nipals", nPcs=2)
\dontshow{stopifnot(sum((fitted(pc) - t(metaboliteData))^2, na.rm=TRUE) < 200)}
}
\author{
  Henning Redestig
}
\references{
  Wold, H. (1966) Estimation of principal components and
  related models by iterative least squares. In
  Multivariate Analysis (Ed., P.R. Krishnaiah), Academic
  Press, NY, 391-420.
}
\seealso{
  \code{prcomp}, \code{princomp}, \code{pca}
}
\keyword{multivariate}

