\name{residuals.pcaRes}
\alias{residuals.pcaRes}
\title{Residuals values from a PCA model.}
\usage{
  \method{residuals}{pcaRes} (object,
    data = completeObs(object), ...)
}
\arguments{
  \item{object}{\code{pcaRes} the \code{pcaRes} object of
  interest.}

  \item{data}{\code{matrix} The data that was used to
  calculate the PCA model (or a different dataset to e.g.
  adress its proximity to the model).}

  \item{...}{Passed on to \code{\link{predict.pcaRes}}.
  E.g. setting the number of used components.}
}
\value{
  A \code{matrix} with the residuals
}
\description{
  This function extracts the residuals values from a pcaRes
  object for the PCA methods SVD, Nipals, PPCA and BPCA
}
\examples{
data(iris)
pcIr <- pca(iris[,1:4])
head(residuals(pcIr, iris[,1:4]))
}
\author{
  Henning Redestig
}
\keyword{multivariate}

