###########################################################################################
##
## bpca <- function(Matrix, nPcs = NULL, completeObs = TRUE, maxSteps = 100, 
##                  verbose = interactive(), ... )
##
## R implementation of a Bayesion PCA missing value estimator.
## After the Matlab script of Shigeyuki OBA (2002  May. 5th)
## See also: http://hawaii.aist-nara.ac.jp/%7Eshige-o/tools/
## Great thanks to them!
##
## Best estimation results can be obtained with the maximum number
## of components (cols(Matrix) - 1). But this is also the computationally
## most expensive option. The complexity is growing with O(n^3) because
## serveral matrix inversions are required. The size of the matrices to invert
## depends on the number of PCs.
## Smaller number of comonents lead also to reasonalble results in
## most cases, but this might depend on your data.
##
## Missing entries in the input data are denoted as 'NA'
## Columns are considered as variables, rows as observations.
##
## Requires:
## Package MASS
##
## Parameters:
## Matrix     - A numeric matrix or data frame. Missing values are
##              denoted as 'NA'
## nPcs       - Number of components used for estimation.
##              The computation time strongly depends on the number of PCs. 
##              The algorithm requires several matrix inversions. The
##              complexity of a matrix inversion is approx. O(n^3), thus
##              complexity increases cubic with the number of used components.
##              The default is (cols(Matrix) - 1).
## completeObs - return completeObs if TRUE
## maxSteps    - Maximum number of estimation steps.
## (optional)    The default is 100
## verbose     - Print some output if TRUE. Default is interactive()
##
##
## Return values:
## pcaRes        - pcaRes is the standart return object for methods
##                 in this package.
##                 Fields are:
##                       pcaRes@completeObs      - estimated complete observations
##                       pcaRes@scores           - estimated scores
##                       pcaRes@loadings         - estimated loadings (eigenvectors)
##                       pcaRes@R2               - the individual R2 values
##                       pcaRes@R2cum            - the cumulative R2 values
##                       pcaRes@sDev             - standart deviation of the scores
##                       pcaRes@nObs             - number of observations (rows of input matrix)
##                       pcaRes@nVar             - number of variables (columns of input matrix)
##                       pcaRes@centered         - boolean, TRUE if centered, FALSE otherwise
##                       pcaRes@center           - row wise mean
##                       pcaRes@varLimit         - NULL, not used here
##                       pcaRes@scaled           - boolean, TRUE if data was scaled, FALSE otherwise
##                       pcaRes@nPcs             - number of principal components (d)
##                       pcaRes@method           - "BPCA"
##                       pcaRes@missing          - Number of NAs in the data
##
##
## Author:    Wolfram Stacklies
##            Max Planck Institut fuer Molekulare Pflanzenphysiologie
##            Golm, Germany
## Date:      04/21/2006
##
## Contact:    wolfram.stacklies@gmail.com
##
###########################################################################################

bpca <- function(Matrix, nPcs = 2, completeObs = TRUE, maxSteps = 100, 
                 verbose = interactive(), ... ) {

    ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix)
    ## And now check if everything is right...
    if ( !checkData(Matrix, verbose = verbose) ) {
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")
    }

    if (nPcs > ncol(Matrix)) {
        stop("more components than matrix columns selected, exiting\n")
    }

    mat <- Matrix

    if (is.null(nPcs)) { nPcs <- ncol(Matrix) - 1 }

    M <- BPCA_initmodel(mat, nPcs)
    tauold <- 1000

    for( step in 1:maxSteps ) {
        M <- BPCA_dostep(M, mat)
        if( step %% 10 == 0 ) {
            tau <- M$tau
            dtau <- abs(log10(tau) - log10(tauold))
            if ( verbose ) {
                cat("Step Number           : ", step, '\n')
                cat("Increase in precision : ", dtau, '\n')
                cat("----------", '\n')
            }
            if (dtau < 1e-4) {
                break
            }
            tauold <- tau
        }
    }
    
    ## Calculate R2cum
    R2cum <- NULL
    centered <- scale(M$yest, center = TRUE, scale = FALSE)
    for (i in 1:nPcs) {
        difference <- centered - ( M$scores[,1:i, drop=FALSE] %*% t(M$PA[,1:i, drop=FALSE]) )
        R2cum <- cbind( R2cum, 1 - ( sum(difference^2) / sum(centered^2) ) )
    }

    ## Calculate R2
    R2 <- vector(length=length(R2cum), mode="numeric")
    R2[1] <- R2cum[1]
    if (nPcs > 1) {
        for (i in 2:nPcs) {
            R2[i] <- R2cum[i] - R2cum[i - 1]
        }
    }

    ####################################################################
    ## Store the values in the pcaRes class, for compatibilty with
    ## other methods provided in the package
    ####################################################################

    result                <- new("pcaRes")

    if (completeObs)
        result@completeObs <- M$yest
    result@center          <- attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center") 
    result@centered        <- TRUE
    result@scaled          <- "none"
    result@scores          <- M$scores 
    result@loadings        <- M$PA
    result@R2cum           <- c(R2cum)
    result@R2              <- R2
    result@sDev            <- apply(M$scores, 2, sd)
    result@nObs            <- nrow(Matrix)
    result@nVar            <- ncol(Matrix)
    result@nPcs            <- nPcs
    result@method          <- "bpca"
    result@missing         <- sum(is.na(Matrix))

    return(result)
}
