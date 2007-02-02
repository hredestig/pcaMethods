##########################################################################################
##
## checkData <- function(data, verbose = FALSE)
##
## Check a given data matrix for consistency with the format
## required for further analysis.
## The data must be a numeric matrix and not contain
## - 'Inf' values
## - 'NaN' values
## - Rows or columns that consist of NA only
##
## Parameters:
##    data - matrix to check
##
## Return values:
##    isValid - TRUE if data passed all tests
##    attributes of isValid:
##        * isNumeric  - TRUE if data is numeric, false otherwise
##        * isInfinite - TRUE if data contains 'Inf' values, false otherwise
##        * isNaN      - TRUE if data contains 'NaN' values, false otherwise
##        * isMatrix   - TRUE if the data is in matrix format, FALSE otherwise
##        * naRows     - TRUE if data contains rows in which all elements are 'NA', 
##                       FALSE otherwise
##        * naCols     - TRUE if data contains columns in which all elements are 'NA', 
##                       FALSE otherwise
##
## Author:  Wolfram Stacklies
##          Max Planck Institute for Molecular Plant Physiology
## Date:    06/28/2006
## Contact: wolfram.stacklies@gmail.com
##
##########################################################################################

checkData <- function(data, verbose = FALSE) {

    ## Initialization
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
            cat("Error: Data is not numeric\n")
    } 

    if ( sum(is.infinite(data) >= 1) ) {
        isInfinite <- TRUE
        isValid <- FALSE
        if (verbose)
            cat("Error: Data contains 'Inf' values\n")
    } 

    if (sum(is.nan(data) >= 1)) {
        isNaN <- TRUE
        isValid <- FALSE
        if (verbose)
            cat("Error: Data contains 'NaN' values\n",
                "Missing values must be denoted by 'NA'.\n")
    } 

    if (!is.matrix(data)) {
        isMatrix <- FALSE
        isValid <- FALSE
        if (verbose)
            cat("Error: data is not a matrix\n",
                "Try to use as.matrix(data)\n")
    }

    ## Check for entire rows that are NA only
    if ( sum(apply(is.na(data), 1, sum) == ncol(data)) >= 1 ) {
        naRows <- TRUE
        isValid <- FALSE
        if (verbose)
            cat("Error: Data contains rows in which all elements are 'NA'\n",
                "Remove them first\n")
    }

    ## Check for entire columns that are NA only
    if ( sum(apply(is.na(data), 2, sum) == nrow(data)) >= 1 ) {
        naCols <- TRUE
        isValid <- FALSE
        if (verbose)
            cat("Error: Data contains columns in which all elements are 'NA'\n",
                "Remove them first\n")
    } 
    
    attr(isValid, "isNumeric")     <- isNumeric
    attr(isValid, "isInfinite")     <- isInfinite
    attr(isValid, "isNaN")         <- isNaN
    attr(isValid, "isMatrix")     <- isMatrix
    attr(isValid, "naRows")     <- naRows
    attr(isValid, "naCols")     <- naCols
    
    return(isValid)
}

