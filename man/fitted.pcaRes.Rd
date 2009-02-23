\name{fitted.pcaRes}
\alias{fitted.pcaRes}
\alias{fitted,pcaRes-method}
\title{Extract fitted values from PCA.}
\description{This function extracts the fitted values from
  a pcaRes object. For PCA methods like SVD, Nipals, PPCA etc
  this is basically just the scores multipled by the loadings, for non-linear
  PCA the original data is propagated through the network to obtain the
  approximated data.}
\usage{fitted.pcaRes(object, data=NULL, nPcs=object@nPcs,...)}
\arguments{
  \item{object}{\code{pcaRes} the \code{pcaRes} object of interest.}
  \item{data}{\code{matrix} For standard PCA methods this can safely be left null to
    get scores x loadings but if set then the scores are obtained by
    projecting provided data onto the loadings. Non-linear PCA is an
    exception, here if data is NULL then data is set to the completeObs
    and propagated through the network.} 
  \item{nPcs}{\code{numeric} The amount of PC's to consider}
  \item{...}{Not passed on anywhere, included for S3 consistency.}
}
\value{A matrix with the fitted values.}
\examples{
data(iris)
pcIr <- pca(iris[,1:4])
head(fitted(pcIr, nPcs=1))
}
\keyword{multivariate}
\author{Henning Redestig <redestig[at]mpimp-golm.mpg.de>}
