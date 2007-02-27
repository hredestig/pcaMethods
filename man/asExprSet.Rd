\name{asExprSet}
\alias{asExprSet}
\title{Convert pcaRes object to an expression set}
\description{
  This function can be used to conveniently replace the 'exprs' slot in
  an \code{exprSet} with the completed data from a \code{pcaRes} object.
}
\usage{
asExprSet(object, exprSet)
}
\arguments{
  \item{object}{\code{pcaRes} -- The object containing the completed
    data.}
  \item{exprSet}{\code{exprSet} -- The object passed on to \code{pca}
    for missing value estimation.}
}
\value{
  An object without missing values of class \code{exprSet}
}
\details{
  This is not a standard \code{as} function as \code{pcaRes} object
  alone not can be converted to an \code{exprSet} (the \code{pcaRes}
  object does not hold any \code{phenoData} for example).
}
\author{Wolfram Stacklies \cr 
  CAS-MPG Partner Institute for Computational Biology, Shanghai, China \cr
  \email{wolfram.stacklies@gmail.com} \cr}
\keyword{multivariate}
