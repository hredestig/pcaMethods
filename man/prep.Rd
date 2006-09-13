\encoding{latin1}
\name{prep}
\alias{prep}
\title{Preprocess a matrix for PCA}
\description{Implements simple preprocessing alternatives for scaling a matrix.}
\usage{prep(object, scale=c("none", "pareto", "vector", "UV"), center=TRUE, ...)}
\arguments{
  \item{object}{Numerical matrix  with (or an object coercible to such) with
    samples in rows and variables as columns. Also takes \code{exprSet} in
    which case the transposed \code{exprs} slot is used.}
  \item{center}{Indicates if the matrix should be mean centred or not.}
  \item{scale}{One of "UV" (unit variance \eqn{a=a/\sigma_{a}})
    "vector" (vector normalisation \eqn{b=b/||b||}), "pareto" or "none"
    to indicate which
    scaling should be used to scale the matrix with \eqn{a} variables
    and \eqn{b} samples.}
  \item{...}{Only used for passing through arguments.}
}
\details{
  Does basically the same as \code{\link{scale}} but adds some
  alternative scaling options. 
}

\value{
  A matrix with attribute "scaled:center" if centring was done.
}
\author{Wolfram Stacklies, Henning Redestig}
\seealso{\code{\link{scale}}}
\examples{
object <- matrix(rnorm(50), nrow=10)
object <- prep(object, scale="vector", center=TRUE)
}
\keyword{manip}

