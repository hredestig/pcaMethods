##' Scaling and centering a matrix.
##'
##' Does basically the same as \code{\link{scale}} but adds some
##' alternative scaling options and functionality for treating
##' pre-processing as part of a model.
##' @title Pre-process a matrix for PCA
##' @param object Numerical matrix (or an object coercible to such)
##' with samples in rows and variables as columns. Also takes
##' \code{ExpressionSet} in which case the transposed expression
##' matrix is used. 
##' @param scale One of "UV" (unit variance \eqn{a=a/\sigma_{a}})
##' "vector" (vector normalisation \eqn{b=b/||b||}), "pareto" (sqrt
##' UV) or "none" to indicate which scaling should be used to scale
##' the matrix with \eqn{a} variables and \eqn{b} samples. Can also be
##' a vector of scales which should be used to scale the
##' matrix. \code{NULL} value is interpreted as \code{"none"}.
##' @param center Either a logical which indicates if the matrix
##' should be mean centred or not, or a vector with averages which
##' should be suntracted from the matrix. \code{NULL} value is
##' interpreted as \code{FALSE}
##' @param eps Minimum variance, variable with lower variance are not
##' scaled and warning is issued instead.
##' @param simple Logical indicating if only the data should be
##' returned or a list with the pre-processing statistics as well.
##' @param reverse Logical indicating  if matrix should be
##' 'post-processed' instead by multiplying each column with its scale
##' and adding the center. In this case, center and scale should be
##' vectors with the statistics (no warning is issued if not, instead
##' output becomes the same as input).
##' @param ... Only used for passing through arguments.
##' @return A pre-processed matrix or a list with
##' \item{center}{a vector with the estimated centers}
##' \item{scale}{a vector with the estimated scales}
##' \item{data}{the pre (or post) processed data}
##' @examples
##' object <- matrix(rnorm(50), nrow=10)
##' res <- prep(object, scale="uv", center=TRUE, simple=FALSE)
##' obj <- prep(object, scale=res$scale, center=res$center)
##' ## same as original
##' sum((object - prep(obj, scale=res$scale, center=res$center, rev=TRUE))^2)
##' @export 
##' @author Henning Redestig
prep <- function(object, scale=c("none", "pareto", "vector", "uv"),
                 center=TRUE, eps=1e-12, simple=TRUE, reverse=FALSE, ...) {
  if(inherits(object, "ExpressionSet"))
    obj <- t(exprs(object))
  else
    obj <- as.matrix(object)

  if(is.null(center))
    center <- FALSE
  if(is.null(scale))
    scale <- "none"

  if(is.logical(center[1])) {
    if(center[1]) 
      center <- colMeans(obj, na.rm=TRUE)
    else
      center <- rep(0, ncol(obj))
  }
  if(length(center) != ncol(obj))
    stop("center do not match matrix dimensions")
  if(!reverse)
    obj <- sweep(obj, 2, center, "-")
  
  if(is.character(scale[1])) {
    scale <- match.arg(scale)
    if(scale == "uv")
      scale <- apply(obj, 2, sd, na.rm=TRUE)
    else if(scale == "none")
      scale <- rep(1, ncol(obj))
    else if(scale == "pareto")
      scale <- sqrt(apply(obj, 2, sd, na.rm=TRUE))
    else if(scale == "vector")
      scale <- apply(obj, 2, function(x) sqrt(sum(x^2, na.rm=TRUE)))
  }
  if(length(scale) != ncol(obj))
    stop("scale vector do not match matrix dimensions")
  
  if (any(scale < eps)) 
    warning(paste("Variance is below eps for", sum(scale < eps),
                  "variables. Not scaling them."))

  scale[scale < eps] <- 1
  if(!reverse)
    obj <- sweep(obj, 2, scale, "/")

  if(reverse) {
    obj <- sweep(obj, 2, scale, "*")
    obj <- sweep(obj, 2, center, "+")
  }
	
  if(inherits(object, "ExpressionSet"))
    exprs(object) <- t(obj)
  else
    object <- obj

  if (simple)
    object
  else
    list(data=object, center=center, scale=scale)
}
