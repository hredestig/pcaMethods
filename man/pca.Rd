\encoding{latin1}
\name{pca}
\alias{pca}
\title{Perform principal component analysis}
\description{Can be used for computing PCA on a numeric
  matrix for visualisation, information extraction and missing value imputation.}
\usage{pca(object, method=c("svd", "nipals", "bpca", "ppca",
"svdImpute", "nlpca", "robustPca"), subset=numeric(),...)}
\arguments{
  \item{object}{Numerical matrix  with (or an object coercible to such) with
    samples in rows and variables as columns. Also takes \code{exprSet} in
    which case the transposed \code{exprs} slot is used.}
  \item{subset}{For convenience one can pass a large matrix but only
    use the variable specified as subset. Can be colnames or indices.}
  \item{method}{One of "svd", "nipals", "bpca", "nlpca" or "ppca".}
  \item{...}{Further arguments to the chosen pca method.}
}
\details{
  This method is wrapper function for the following set of pca methods:
  \describe{
    \item{svd:}{Uses classical \code{prcomp}. See documentation for \code{\link{svdPca}}.}
    \item{nipals:}{An iterative method capable of handling small amounts
      of missing values. See documentation for \code{\link{nipalsPca}}.}
    \item{bpca:}{An iterative method using a Bayesian model to handle missing values. See documentation for \code{\link{bpca}}.}
    \item{ppca:}{An iterative method using a probabilistic model to handle missing values. See documentation for \code{\link{ppca}}.}
    \item{svdImpute:}{Uses expectation maximation to perform SVD PCA on incomplete data. See documentation for \code{\link{svdImpute}}.}
  }
  Extra arguments usually given to this function include:
  \describe{
    \item{nPcs:}{The amount of principal components to extract}
  }
}

\value{
  A \code{pcaRes} object. Or a list containing a pcaRes object as first and an
  exprSet object as second entry if the input was of type exprSet.
}
\references{
  Wold, H. (1966) Estimation of principal components and related models by
  iterative least squares. In Multivariate Analysis (Ed.,
  P.R. Krishnaiah), Academic Press, NY, 391-420.
  
  Shigeyuki Oba, Masa-aki Sato, Ichiro Takemasa,
  Morito Monden, Ken-ichi Matsubara and Shin Ishii.
  A Bayesian missing value estimation method for gene expression profile
  data. \emph{Bioinformatics, 19(16):2088-2096, Nov 2003}.

  Troyanskaya O. and Cantor M. and Sherlock G. and Brown
  P. and Hastie T. and Tibshirani R. and Botstein D. and Altman RB.
  - Missing value estimation methods for DNA microarrays.
  \emph{Bioinformatics. 2001 Jun;17(6):520-5}.
}
\author{Wolfram Stacklies, Henning Redestig}
\seealso{\code{\link{prcomp}}, \code{\link{princomp}},
  \code{\link{nipalsPca}}, \code{\link{svdPca}}}
\examples{
data(iris)
## Usually some kind of scaling is appropriate
pcIr <- pca(iris[,1:4], nPcs = 2, method="nipals")
pcIr <- pca(iris[,1:4], nPcs = 2, method="svd")
## Get a short summary on the calculated model
summary(pcIr)
## Scores and loadings plot
slplot(pcIr, sl=as.character(iris[,5]))
}
\keyword{multivariate}
