\name{BPCA_dostep}
\alias{BPCA_dostep}
\title{Do BPCA estimation step}
\usage{
  BPCA_dostep(M, y)
}
\arguments{
  \item{M}{Data structure containing all needed
  information. See the source documentation of
  BPCA_initmodel for details}

  \item{y}{Numeric original data matrix}
}
\value{
  Updated version of the data structure
}
\description{
  The function contains the actual implementation of the
  BPCA component estimation. It performs one step of the
  BPCA EM algorithm.  It is called 'maxStep' times from
  within the main loop in BPCAestimate.
}
\details{
  This function is NOT intended to be run standalone.
}
\author{
  Wolfram Stacklies
}

