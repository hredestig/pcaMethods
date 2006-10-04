\name{svdPca}
\alias{svdPca}
\title{Perform principal component analysis using singular value decomposition}
\description{A wrapper function for R's standard function \code{prcomp}.
  Delivers the result as a \code{pcaRes} method for compatibility with
  the rest of the pcaMethods package.}
\usage{svdPca(Matrix, nPcs=2, center=TRUE, completeObs=FALSE, varLimit=1,...)}
\arguments{
  \item{Matrix}{Numerical matrix samples in rows and variables as columns.}
  \item{nPcs}{Number of components that should be extracted.}
  \item{center}{Center the data column wise if TRUE}
  \item{completeObs}{Return the complete observations. This exisits for compatibility
  only, as svdPca cannot missing values. If set TRUE the input matrix will be returned
  in the \code{completeObs} field.}
  \item{varLimit}{Optionally the ratio of variance that should be
    explained. \code{nPcs} is ignored if varLimit < 1}
  \item{...}{Only used for passing through arguments.}
}
\details{
  svdPca can preferrably be called using \code{pca(object, method="svd")}.
}

\value{
  A \code{pcaRes} object.
}
\author{Henning Redestig}
\seealso{\code{prcomp}, \code{princomp}, \code{pca}}
\examples{
data(iris)
pcIr <- svdPca(iris[,1:4], nPcs=2)
}
\keyword{multivariate}
