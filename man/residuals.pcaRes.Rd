\name{residuals.pcaRes}
\alias{residuals.pcaRes}
\alias{residuals,pcaRes-method}
\title{Residuals values from a PCA model.}
\description{This function extracts the residuals values from
  a pcaRes object for the PCA methods SVD, Nipals, PPCA and BPCA
}

\usage{residuals.pcaRes(object, data, nPcs=object@nPcs, ...)}
\arguments{
  \item{object}{\code{pcaRes} the \code{pcaRes} object of interest.}
  \item{data}{\code{matrix} The data that was used to calculate the PCA
    model (or a different dataset to e.g. adress its proximity to the model).} 
  \item{nPcs}{\code{numeric} The amount of PC's to consider}
  \item{...}{Not passed on anywhere, included for S3 consistency.}
}
\value{A \code{matrix} with the residuals
}
\keyword{multivariate}
\examples{
data(iris)
pcIr <- pca(iris[,1:4])
head(residuals(pcIr, iris[,1:4]))
}
\author{Henning Redestig <redestig[at]psc.riken.jp>}
