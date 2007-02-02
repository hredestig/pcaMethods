\encoding{latin1}
\name{Q2}
\alias{Q2}
\title{Perform internal cross-validation for PCA}
\description{Internal cross-validation can be used for estimating the
  level of structure in a data set and to optimise the choice of number
  of principal components.}
\usage{Q2(object, originalData, nPcs=object@nPcs, fold=5, nruncv=10,
segments=NULL, verbose=interactive())}
\arguments{
\item{object}{A \code{pcaRes} object (result from previous PCA
analysis.)}
\item{originalData}{The matrix used to obtain the pcaRes object}
\item{nPcs}{The amount of principal components to estimate Q2 for.}
\item{fold}{The amount of groups to divide the data in.}
\item{nruncv}{The amount of times to repeat the whole cross-validation}
\item{segments}{\code{list} A predefined list where each element is the
  set of indices to leave out. Note that if this is provided, Q2 becomes
  deterministic (if the PCA is deterministic of course).}
\item{verbose}{\code{boolean} If TRUE Q2 outputs a primitive progress bar.}
}
\details{
This method calculates \eqn{Q^2} for a PCA model. This is the predictory
version of \eqn{R^2} and can be interpreted as the ratio of
variance in a left out data chunk that can be estimated by the PCA
model. Poor (low) \eqn{Q^2} means that the PCA model only describes noise
and that the model is unrelated to the true data structure. The
definition of \eqn{Q^2} is:

\deqn{Q^2 = 1 - \frac{\sum_{i}^{k}\sum_{j}^{n}(x -
\hat{x})^2}{\sum_{i}^{k}\sum_{j}^{n}x^2}}{Q^2 = 1 - sum_i^k sum_j^n (x -
\hat{x})^2 / \sum_i^k \sum_j^n(x^2)}

for the matrix \eqn{x} which has \eqn{n} rows and \eqn{k}
columns. For a given amount of PC's x is estimated as \eqn{\hat{x} =
TP'} (T are scores and P are loadings). Though this defines the
leave-one-out cross-validation this is 
not what is performed if fold is less than the amount of rows and/or
columns.

Diagonal rows of elements in the matrix are deleted and the
re-estimated. You can choose your own segmentation as well make sure no
complete row or column is lost.
}

\value{
  A matrix with \eqn{Q^2} estimates.
}
\references{
  Wold, H. (1966) Estimation of principal components and related models by
  iterative least squares. In Multivariate Analysis (Ed.,
  P.R. Krishnaiah), Academic Press, NY, 391-420.
}
\author{Wolfram Stacklies, Henning Redestig}
\seealso{\code{\link{pca}}}
\examples{
data(iris)
pcIr <- pca(iris[,1:4], nPcs=2, scale="UV", method="ppca")
#can only get Q2 estimats for the two first PC's
q2 <- Q2(pcIr, iris[,1:4], nruncv=2)
#Typically Q2 increases only very slowly after the optimal amount of PC's
boxplot(q2~row(q2), xlab="Amount of PC's", ylab=expression(Q^2))
}
\keyword{multivariate}
