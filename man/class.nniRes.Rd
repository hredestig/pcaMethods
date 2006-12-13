\name{nniRes}
\docType{class}
\alias{nniRes}
\alias{nniRes-class}
\alias{print,nniRes-method}

\title{Class for representing a nearest neighbour imputation result}

\description{This is a class representation of nearest neighbour imputation (nni)
  result}

\section{Creating Objects}{
  \code{new("nniRes",
    completeObs=[the estimated complete observations],
    k=[cluster size], nObs=[amount of observations], nVar=[amount of variables],
    scaled=[which method was used to scale data], centered=[was data centered],
    method=[method used to perform clustering], subset=[subset of variables of data used],
    missing=[amount of NAs], center=[original means])}
}

\section{Slots}{
  \describe{
    \item{completeObs}{"matrix", the estimated complete observations}
    \item{nObs}{"numeric", amount of observations}
    \item{nVar}{"numeric", amount of variables}
    \item{centered}{"logical", data was centered or not}
    \item{center}{"numeric", the original variable centers}
    \item{subset}{"numeric", the subset of variables used}
    \item{scaled}{"character", the scaled method  used}
    \item{k}{"numeric", cluster size}
    \item{method}{"character", the method used to perform the clustering}
    \item{missing}{"numeric", the total amount of missing values in
      original data}
  }
}
\section{Methods}{
  \describe{
    \item{print}{Print function}
  }
}
\keyword{classes}
