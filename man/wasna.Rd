\name{wasna}
\alias{wasna}
\alias{wasna,pcaRes-method}
\title{Get a matrix with indicating the elements that were missing in the
input data. Convenient for estimating imputation performance.}
\usage{
  wasna(object, ...)
}
\arguments{
  \item{object}{pcaRes object}

  \item{...}{Not used}
}
\value{
  A matrix with logicals
}
\description{
  Get a matrix with indicating the elements that were
  missing in the input data. Convenient for estimating
  imputation performance.
}
\examples{
data(metaboliteData)
data(metaboliteDataComplete)
result <- pca(metaboliteData, nPcs=2)
plot(completeObs(result)[wasna(result)], metaboliteDataComplete[wasna(result)])
}
\author{
  Henning Redestig
}

