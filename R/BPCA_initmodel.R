######################################################################################
##
## Model initialization for Bayesian PCA.
## This function is NOT inteded to be run separately!
##
## The function calculates the initial Eigenvectors by use
## of SVD from the complete rows.
## The data structure M is created and initial values are 
## assigned.
##
## Parameters:
## y        - numeric matrix containing missing values.
##          Missing values are denoted as 'NA'
## q        - Number of components used for estimation
##
## Return values:
## M             - Data Structure containing all runtime parameters.
##                 These are:
##                 rows          : Row number of input matrix
##                 cols          : Column number of input matrix
##                 comps         : Number of components to use
##                 yest          : (working variable) current estimate of complete data
##                 row_miss      : (Array) Indizes of rows containing missing values
##                 row_nomiss    : (Array) Indices of complete rows (such with no missing values)
##                 nans          : Matrix of same size as input data.
##                                 TRUE if input == NA, false otherwise
##                 mean          : Column wise data mean
##                 PA            : (d x k) Estimated principal axes (eigenvectors, loadings)
##                                 The matrix ROWS are the vectors
##                 tau           : Estimated precision of the residual error
##          scores    : Estimated scores
##                 
##                 Further elements are:
##                 galpha0, balpha0, alpha, gmu0, btau0, gtau0, SigW
##                 These are working variables or constants.
##
##
## Author:       Wolfram Stacklies
##               Max Planck Institut fuer Molekulare Pflanzenphysiologie
##               Golm, Germany
## Date:         04/21/2006
##
## Contact:      wolfram.stacklies@gmail.com
##
######################################################################################

BPCA_initmodel <- function(y, components) {
    ## Initialization, write static parameters to the central
    M <- NULL ## data structure
    M$rows <- nrow(y) ## Row number
    M$cols <- ncol(y) ## Number of components to use for estimation 
    M$comps <- components ## Column number
    M$yest <- y ## Original data, NAs are set to 0 later on

    ## Find rows with missing values, etc...
    M$nans <- is.na(y)
    temp <- apply(M$nans, 1, sum)
    M$row_nomiss <- which(temp == 0)
    M$row_miss <- which(temp != 0)
    M$yest[M$nans] <- 0
    M$scores <- NULL

    ## Get the SVD of the complete rows
    covy <- cov(M$yest)
    values <- svd(covy, components, components)
    U <- values[[2]]
    S <- diag( values[[1]][1:components], nrow = components, ncol = components)
    V <- values[[3]]

    ## M$mean: column wise mean of the original data
    M$mean <- matrix(0, 1, M$cols)
    for(j in 1:M$cols) {
        idx <- which(!is.na(y[,j]))
        M$mean[j] <- mean(y[idx,j])
    }

    M$PA <- U %*% sqrt(S)
    M$tau <- 1 / ( sum(diag(covy)) - sum(diag(S)) )
    
    ## Constants etc
    taumax <- 1e10
    taumin <- 1e-10
    M$tau <- max( min(M$tau, taumax), taumin )

    M$galpha0 <- 1e-10
    M$balpha0 <- 1
    M$alpha <- (2 * M$galpha0 + M$cols) / (M$tau * diag(t(M$PA) %*% M$PA) + 2 * M$galpha0 / M$balpha0)

    M$gmu0 <- 0.001

    M$btau0 <- 1
    M$gtau0 <- 1e-10
    M$SigW <- diag(components)
    return(M)
}
