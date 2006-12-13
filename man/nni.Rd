\encoding{latin1}
\name{nni}
\alias{nni}
\title{Nearest neighbour imputation}
\description{Wrapper function for imputation methods based on nearest
neighbour clustering. Currently llsImpute only.}
\usage{nni(object, method=c("llsImpute"), subset=numeric(),...)}
\arguments{
  \item{object}{Numerical matrix  with (or an object coercible to such) with
    samples in rows and variables as columns. Also takes \code{exprSet} in
    which case the transposed \code{exprs} slot is used.}
  \item{subset}{For convenience one can pass a large matrix but only
    use the variable specified as subset. Can be colnames or indices.}
  \item{method}{Currently "llsImpute" only.}
  \item{...}{Further arguments to the chosen method.}
}
\details{
  This method is wrapper function to llsImpute, See documentation for \code{link{llsImpute}}
  Extra arguments usually given to this function include:
}

\value{
  A \code{clusterRes} object. Or a list containing a clusterRes object as first and an
  exprSet object as second entry if the input was of type exprSet.
}
\author{Wolfram Stacklies}
\seealso{ \code{\link{llsImpute}}, \code{\link{pca}}}
\examples{
data(metaboliteData)
llsRes <- nni(metaboliteData, k=6, method="llsImpute", allGenes=TRUE)
}
\keyword{multivariate}
