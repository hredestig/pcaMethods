\name{nipalsPca}
\alias{nipalsPca}
\title{Perform principal component analysis using the Non-linear iterative
  partial least squares (NIPALS) algorithm.}
\description{Can be used for computing PCA on a numeric
  matrix using either the NIPALS
  algorithm which is an iterative approach for estimating the principal
  components extracting them one at a time. NIPALS can handle a small
  amount of missing values.}
\usage{nipalsPca(object, nPcs=2, varLimit=1, maxSteps=5000, limit=1e-6, verbose=interactive(),...)}
\arguments{
  \item{object}{Numerical matrix samples in rows and variables as
    columns.}
  \item{maxSteps}{Defines how many iterations can be done before the
    algorithm should abort (happens almost exclusively when there were
    some wrong in the input data).}
  \item{limit}{The limit condition for judging if the algorithm has
    converged or not, specifically if a new iteration is done if
    \eqn{(T_{old} - T)^T(T_{old} - T) > \code{limit}}.}
  \item{nPcs}{Number of components that should be extracted.}
  \item{varLimit}{Optionally the ratio of variance that should be
    explained. \code{nPcs} is ignored if varLimit < 1}
  \item{verbose}{Show simple progress information.}
  \item{...}{Only used for passing through arguments.}
}
\details{
  This method is incredibly slow and should not be used on larger
  matrices. The power in missing value imputation is also quite
  disputable.
}

\value{
  A \code{pcaRes} object.
}
\references{
  Wold, H. (1966) Estimation of principal components and related models by
  iterative least squares. In Multivariate Analysis (Ed.,
  P.R. Krishnaiah), Academic Press, NY, 391-420.
}
\author{Henning Redestig}
\seealso{\code{prcomp}, \code{princomp}, \code{pca}}
\examples{
data(iris)
pcIr <- nipalsPca(iris[,1:4], nPcs=2)
}
\keyword{multivariate}
