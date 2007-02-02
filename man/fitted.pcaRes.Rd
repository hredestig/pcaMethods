\name{fitted.pcaRes}
\alias{fitted.pcaRes}
\alias{fitted,pcaRes-method}
\title{Extract fitted values from PCA.}
\description{This function extracts the fitted values from
  apcaResobject. For PCA methods like SVD, Nipals, PPCA etc
  thisisbasically just the scores multipled by the loadings,forNon-linear
  PCA the original data is propagated throughthenetwork to obtain the
  approximated data.}
\usage{fitted.pcaRes(object, data=NULL, nPcs=object@nPcs)}
\arguments{
  \item{object}{the \code{pcaRes} object of interest.}
  \item{data}{For standard PCA methods this can safely be left null to
    get scores x loadings but if set then the scores are obtained by
    projecting provided data onto the loadings. Non-linear PCA is an
    exception, here if data is NULL then data is set to the completeObs
    and propagated through the network.} 
  \item{nPcs}{The amount of PC's to consider}
}
\value{A matrix with the fitted values.}
\keyword{multivariate}
\author{Henning Redestig <redestig[at]mpimp-golm.mpg.de>}
