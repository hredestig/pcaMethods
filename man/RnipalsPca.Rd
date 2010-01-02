\name{RnipalsPca}
\alias{RnipalsPca}
\title{NIPALS PCA implemented in R}
\usage{RnipalsPca(Matrix, nPcs=2, center=TRUE, completeObs=TRUE, varLimit=1,
    maxSteps=5000, threshold=1e-06, verbose=interactive(), ...)}
\description{PCA by non-linear iterative partial least squares}
\details{Can be used for computing PCA on a numeric matrix using either the
NIPALS algorithm which is an iterative approach for estimating the
principal components extracting them one at a time. NIPALS can
handle a small amount of missing values. It is not recommended to
use this function directely but rather to use the pca() wrapper
function. There is a C++ implementation given as \code{nipalsPca} which is faster.}
\value{A \code{pcaRes} object.}
\references{Wold, H. (1966) Estimation of principal components and related models by
iterative least squares. In Multivariate Analysis (Ed.,
P.R. Krishnaiah), Academic Press, NY, 391-420.}
\author{Henning Redestig}
\seealso{\code{prcomp}, \code{princomp}, \code{pca}}
\keyword{multivariate}
\arguments{\item{Matrix}{Numerical matrix samples in rows and variables as columns.}
\item{nPcs}{Number of components that should be extracted.}
\item{center}{Mean center the data column wise if set TRUE}
\item{completeObs}{Return the estimated complete observations. This is
the input Matrix with NA values replaced by the estimated values.}
\item{varLimit}{Optionally the ratio of variance that should be
explained. \code{nPcs} is ignored if varLimit < 1}
\item{maxSteps}{Defines how many iterations can be done before
algorithm should abort (happens almost exclusively when there were
some wrong in the input data).}
\item{threshold}{The limit condition for judging if the algorithm has
converged or not, specifically if a new iteration is done if
\eqn{(T_{old} - T)^T(T_{old} - T) > \code{limit}}.}
\item{verbose}{Show simple progress information.}
\item{...}{Only used for passing through arguments.}}
\examples{data(iris)
mat <- as.matrix(iris[,1:4])
system.time(pcIr <- RnipalsPca(mat, method="rnipals", nPcs=3))
system.time(pcIr <- nipalsPca(mat, nPcs=3))}
