\name{checkData}
\alias{checkData}
\title{Do some basic checks on a given data matrix}
\usage{
  checkData(data, verbose = FALSE)
}
\arguments{
  \item{data}{\code{matrix} -- Data to check.}

  \item{verbose}{\code{boolean} -- If TRUE, the function
  prints messages whenever an error in the data set is
  found.}
}
\value{
  \item{isValid}{\code{boolean} -- TRUE if no errors were
  found, FALSE otherwise.  isValid contains a set of
  attributes, these are: \itemize{ \item isNumeric - TRUE
  if data is numeric, false otherwise \item isInfinite -
  TRUE if data contains 'Inf' values, false otherwise \item
  isNaN - TRUE if data contains 'NaN' values, false
  otherwise \item isMatrix - TRUE if the data is in matrix
  format, FALSE otherwise \item naRows - TRUE if data
  contains rows in which all elements are 'NA', FALSE
  otherwise \item naCols - TRUE if data contains columns in
  which all elements are 'NA', FALSE otherwise }}
}
\description{
  Check a given data matrix for consistency with the format
  required for further analysis. The data must be a numeric
  matrix and not contain: \itemize{ \item Inf values \item
  NaN values \item Rows or columns that consist of NA only
  }
}
\author{
  Wolfram Stacklies
}
\keyword{multivariate}

