\name{asExprSet}
\alias{asExprSet}
\title{Convert pcaRes object to an expression set}
\usage{asExprSet(object, exprSet)}
\description{This function can be used to conveniently replace the expression
matrix in an \code{ExpressionSet} with the completed data from a
\code{pcaRes} object.}
\details{This is not a standard \code{as} function as \code{pcaRes}
object alone not can be converted to an \code{ExpressionSet} (the
\code{pcaRes} object does not hold any \code{phenoData} for
example).}
\value{An object without missing values of class \code{ExpressionSet}.}
\author{Wolfram Stacklies \cr CAS-MPG Partner Institute for
Computational Biology, Shanghai, China}
\keyword{multivariate}
\arguments{\item{object}{\code{pcaRes} -- The object containing the completed
data.}
\item{exprSet}{\code{ExpressionSet} -- The object passed on to
\code{pca} for missing value estimation.}}
