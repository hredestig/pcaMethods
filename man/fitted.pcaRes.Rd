\name{fitted.pcaRes}
\alias{fitted.pcaRes}
\title{Extract fitted values from PCA.}
\usage{
  \method{fitted}{pcaRes} (object, data = NULL,
    nPcs = nP(object), pre = TRUE, post = TRUE, ...)
}
\arguments{
  \item{object}{the \code{pcaRes} object of interest.}

  \item{data}{For standard PCA methods this can safely be
  left null to get scores x loadings but if set, then the
  scores are obtained by projecting provided data onto the
  loadings.  If data contains missing values the result
  will be all NA. Non-linear PCA is an exception, here if
  data is NULL then data is set to the completeObs and
  propaged through the network.}

  \item{nPcs}{The number of PC's to consider}

  \item{pre}{pre-process \code{data} based on the
  pre-processing chosen for the PCA model}

  \item{post}{unpre-process the final data (add the center
  back etc to get the final estimate)}

  \item{...}{Not used}
}
\value{
  A matrix representing the fitted data
}
\description{
  Fitted values of a PCA model
}
\details{
  This function extracts the fitted values from a
  pcaResobject. For PCA methods like SVD, Nipals, PPCA etc
  this is basically just the scores multipled by the
  loadings and adjusted for pre-processing. for non-linear
  PCA the original data is propagated through the network
  to obtain the approximated data.
}
\examples{
pc <- pca(iris[,1:4], nPcs=4, center=TRUE, scale="uv")
sum( (fitted(pc) - iris[,1:4])^2 )
}
\author{
  Henning Redestig
}
\keyword{multivariate}

