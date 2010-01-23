\name{nni}
\alias{nni}
\title{Nearest neighbour imputation}
\usage{nni(object, method=c("llsImpute"), subset=numeric(), ...)}
\description{Wrapper function for imputation methods based on nearest neighbour
clustering. Currently llsImpute only.}
\details{This method is wrapper function to llsImpute, See documentation
for \code{link{llsImpute}}.}
\value{A \code{clusterRes} object. Or a list containing a
clusterRes object as first and an ExpressionSet object as second
entry if the input was of type ExpressionSet.}
\seealso{\code{\link{llsImpute}}, \code{\link{pca}}}
\keyword{multivariate}
\author{Wolfram Stacklies}
\arguments{\item{object}{Numerical matrix with (or an object coercible to
such) with samples in rows and variables as columns. Also takes
\code{ExpressionSet} in which case the transposed expression
matrix is used.}
\item{method}{For convenience one can pass a large matrix but only
use the variable specified as subset. Can be colnames or indices.}
\item{subset}{Currently "llsImpute" only.}
\item{...}{Further arguments to the chosen method.}}
\examples{data(metaboliteData)
llsRes <- nni(metaboliteData, k=6, method="llsImpute", allGenes=TRUE)}
