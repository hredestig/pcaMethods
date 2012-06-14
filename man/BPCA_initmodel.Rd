\name{BPCA_initmodel}
\alias{BPCA_initmodel}
\title{Initialize BPCA model}
\usage{
  BPCA_initmodel(y, components)
}
\arguments{
  \item{y}{numeric matrix containing missing values.
  Missing values are denoted as 'NA'}

  \item{components}{Number of components used for
  estimation}
}
\value{
  List containing \item{rows}{Row number of input matrix}
  \item{cols}{Column number of input matrix}
  \item{comps}{Number of components to use}
  \item{yest}{(working variable) current estimate of
  complete data} \item{row_miss}{(Array) Indizes of rows
  containing missing values} \item{row_nomiss}{(Array)
  Indices of complete rows (such with no missing values)}
  \item{nans}{Matrix of same size as input data. TRUE if
  \code{input == NA}, false otherwise} \item{mean}{Column
  wise data mean} \item{PA}{ (d x k) Estimated principal
  axes (eigenvectors, loadings) The matrix ROWS are the
  vectors} \item{tau}{Estimated precision of the residual
  error} \item{scores}{ Estimated scores} Further elements
  are: galpha0, balpha0, alpha, gmu0, btau0, gtau0, SigW.
  These are working variables or constants.
}
\description{
  Model initialization for Bayesian PCA. This function is
  NOT inteded to be run separately!
}
\details{
  The function calculates the initial Eigenvectors by use
  of SVD from the complete rows.  The data structure M is
  created and initial values are assigned.
}
\author{
  Wolfram Stacklies
}

