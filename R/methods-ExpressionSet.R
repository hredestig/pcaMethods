##' This function can be used to conveniently replace the expression
##' matrix in an \code{ExpressionSet} with the completed data from a
##' \code{pcaRes} object.
##'
##' This is not a standard \code{as} function as \code{pcaRes}
##' object alone not can be converted to an \code{ExpressionSet} (the
##' \code{pcaRes} object does not hold any \code{phenoData} for
##' example).
##' @title Convert pcaRes object to an expression set
##' @param object \code{pcaRes} -- The object containing the completed
##' data.
##' @param exprSet \code{ExpressionSet} -- The object passed on to
##' \code{pca} for missing value estimation.  
##' @return An object without missing values of class \code{ExpressionSet}.
##' @export 
##' @author Wolfram Stacklies \cr CAS-MPG Partner Institute for
##' Computational Biology, Shanghai, China
##' @keywords multivariate
asExprSet <- function(object, exprSet) {
  if(!inherits(exprSet, "ExpressionSet"))
    stop("Parameter exprSet must be of type ExpressionSet")
  if(!inherits(object, "pcaRes") & !inherits(object, "nniRes"))
    stop("Parameter object must be either of type pcaRes or nniRes")
  if (is.null(completeObs(object)))
    stop("completeObs(object) is NULL, exiting")
  if(!all(dim(exprs(exprSet)) == dim(t(completeObs(object)))))
    stop("Dimensions of exprs(exprSet) and completeObs(object) do not match. 
Did you really do missing value estimation using this ExpressionSet object?")
  exprs(exprSet) <- t(completeObs(object))
  return(exprSet)
}
