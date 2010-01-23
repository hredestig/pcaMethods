\name{wasna}
\title{Get a matrix with indicating the elements that were missing in the...}
\usage{wasna(object, ...)}
\description{Get a matrix with indicating the elements that were missing in the
input data. Convenient for estimating imputation performance.}
\alias{wasna}
\alias{wasna,pcaRes-method}
\value{A matrix with logicals}
\author{Henning Redestig}
\arguments{\item{object}{pcaRes object}
\item{...}{Not used}}
\examples{data(metaboliteData)
data(metaboliteDataComplete)
result <- pca(metaboliteData, nPcs=2)
plot(completeObs(result)[wasna(result)], metaboliteDataComplete[wasna(result)])}
