\docType{class}
\name{nniRes}
\alias{nniRes}
\alias{nniRes-class}
\title{Class for representing a nearest neighbour imputation result}
\description{
  This is a class representation of nearest neighbour
  imputation (nni) result
}
\details{
  \bold{Creating Objects}\cr \code{new("nniRes",
  completeObs=[the estimated complete observations],
  k=[cluster size], nObs=[amount of observations],
  nVar=[amount of variables], centered=[was the data
  centered befor running LLSimpute], center=[original
  means], method=[method used to perform clustering],
  missing=[amount of NAs])}

  \bold{Slots}\cr \describe{ \item{completeObs}{"matrix",
  the estimated complete observations}
  \item{nObs}{"numeric", amount of observations}
  \item{nVar}{"numeric", amount of variables}
  \item{correlation}{"character", the correlation method
  used (pearson, kendall or spearman)}
  \item{centered}{"logical", data was centered or not}
  \item{center}{"numeric", the original variable centers}
  \item{k}{"numeric", cluster size}
  \item{method}{"character", the method used to perform the
  clustering} \item{missing}{"numeric", the total amount of
  missing values in original data} }

  \bold{Methods}\cr \describe{ \item{print}{Print function}
  }
}
\author{
  Wolfram Stacklies
}
\keyword{classes}

