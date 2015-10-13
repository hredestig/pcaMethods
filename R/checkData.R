##' Check a given data matrix for consistency with the format
##' required for further analysis.
##' The data must be a numeric matrix and not contain:
##' \itemize{
##' \item Inf values
##' \item NaN values
##' \item Rows or columns that consist of NA only
##' }
##' @title Do some basic checks on a given data matrix
##' @param data \code{matrix} -- Data to check.
##' @param verbose  \code{boolean} -- If TRUE, the function prints
##' messages whenever an error in the data set is found.
##' @return \item{isValid}{\code{boolean} -- TRUE if no errors were
##' found, FALSE otherwise.  isValid contains a set of attributes,
##' these are: \itemize{ \item isNumeric - TRUE if data is numeric,
##' false otherwise \item isInfinite - TRUE if data contains 'Inf'
##' values, false otherwise \item isNaN - TRUE if data contains 'NaN'
##' values, false otherwise \item isMatrix - TRUE if the data is in
##' matrix format, FALSE otherwise \item naRows - TRUE if data
##' contains rows in which all elements are 'NA',  FALSE otherwise
##' \item naCols - TRUE if data contains columns in which all elements
##' are 'NA',  FALSE otherwise }}
##' @keywords multivariate
##' @export 
##' @author Wolfram Stacklies
checkData <- function(data, verbose = FALSE) {
  isValid <- TRUE

  isNumeric <- TRUE
  isInfinite <- FALSE
  isNaN <- FALSE
  isMatrix <- TRUE
  naRows <- FALSE
  naCols <- FALSE

  if (!is.numeric(data)) {
    isNumeric <- FALSE
    isValid <- FALSE
    if (verbose)
      message("Error: Data is not numeric")
  } 

  if ( sum(is.infinite(data) >= 1) ) {
    isInfinite <- TRUE
    isValid <- FALSE
    if (verbose)
      message("Error: Data contains 'Inf' values")
  } 

  if (sum(is.nan(data) >= 1)) {
    isNaN <- TRUE
    isValid <- FALSE
    if (verbose)
      message("Error: Data contains 'NaN' values. Missing values must be denoted by 'NA'")
  } 

  if (!is.matrix(data)) {
    isMatrix <- FALSE
    isValid <- FALSE
    if (verbose)
      message("Error: data is not a matrix. Try to use as.matrix(data)")
  }

  ## Check for entire rows that are NA only
  if (sum(apply(is.na(data), 1, sum) == ncol(data)) >= 1 ) {
    naRows <- TRUE
    isValid <- FALSE
    if (verbose)
      message("Error: Data contains rows in which all elements are 'NA'. Remove them first")
  }

  ## Check for entire columns that are NA only
  if (sum(apply(is.na(data), 2, sum) == nrow(data)) >= 1 ) {
    naCols <- TRUE
    isValid <- FALSE
    if (verbose)
      message("Error: Data contains columns in which all elements are 'NA'. Remove them first")
  } 
  
  attr(isValid, "isNumeric") <- isNumeric
  attr(isValid, "isInfinite") <- isInfinite
  attr(isValid, "isNaN") <- isNaN
  attr(isValid, "isMatrix") <- isMatrix
  attr(isValid, "naRows") <- naRows
  attr(isValid, "naCols") <- naCols
  return(isValid)
}

