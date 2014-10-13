\name{prep}
\alias{prep}
\title{Pre-process a matrix for PCA}
\usage{
  prep(object, scale = c("none", "pareto", "vector", "uv"),
    center = TRUE, eps = 1e-12, simple = TRUE,
    reverse = FALSE, ...)
}
\arguments{
  \item{object}{Numerical matrix (or an object coercible to
  such) with samples in rows and variables as columns. Also
  takes \code{ExpressionSet} in which case the transposed
  expression matrix is used.}

  \item{scale}{One of "UV" (unit variance
  \eqn{a=a/\sigma_{a}}) "vector" (vector normalisation
  \eqn{b=b/||b||}), "pareto" (sqrt UV) or "none" to
  indicate which scaling should be used to scale the matrix
  with \eqn{a} variables and \eqn{b} samples. Can also be a
  vector of scales which should be used to scale the
  matrix. \code{NULL} value is interpreted as
  \code{"none"}.}

  \item{center}{Either a logical which indicates if the
  matrix should be mean centred or not, or a vector with
  averages which should be suntracted from the matrix.
  \code{NULL} value is interpreted as \code{FALSE}}

  \item{eps}{Minimum variance, variable with lower variance
  are not scaled and warning is issued instead.}

  \item{simple}{Logical indicating if only the data should
  be returned or a list with the pre-processing statistics
  as well.}

  \item{reverse}{Logical indicating if matrix should be
  'post-processed' instead by multiplying each column with
  its scale and adding the center. In this case, center and
  scale should be vectors with the statistics (no warning
  is issued if not, instead output becomes the same as
  input).}

  \item{...}{Only used for passing through arguments.}
}
\value{
  A pre-processed matrix or a list with \item{center}{a
  vector with the estimated centers} \item{scale}{a vector
  with the estimated scales} \item{data}{the pre (or post)
  processed data}
}
\description{
  Scaling and centering a matrix.
}
\details{
  Does basically the same as \code{\link{scale}} but adds
  some alternative scaling options and functionality for
  treating pre-processing as part of a model.
}
\examples{
object <- matrix(rnorm(50), nrow=10)
res <- prep(object, scale="uv", center=TRUE, simple=FALSE)
obj <- prep(object, scale=res$scale, center=res$center)
## same as original
sum((object - prep(obj, scale=res$scale, center=res$center, rev=TRUE))^2)
}
\author{
  Henning Redestig
}

