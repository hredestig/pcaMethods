##########################################################################################
##
## ppca <- function(Matrix, nPcs, center = TRUE, completeObs = TRUE, seed = NA, ... )
##
## Implements probabilistic PCA for data with missing values,
## using a factorizing distribution over hidden states and hidden
## observations.
##
## This script is a port from matlab to R. The original matlab
## version was contributed by J.J Verbeek in 2002, see also
## http://www.science.uva.nl/~jverbeek
## Thanks a lot!
##
## Missing entries in the input data are denoted as 'NA'
## Columns are considered as variables, rows as observations.
##
## Parameters:
## Matrix      - A numeric matrix or data frame. Missing values are denoted
##               as 'NA'
## nPcs        - Number of principal components to calculate
## center      - Mean center the data, if TRUE
## completeObs - Return estimated complete observations if TRUE
## seed        - Set the seed for the random number generator. PPCA creates
##               fills the initial loading matrix with random numbers chosen from a normal
##               distribution. Thus results may vary slightly. Set the seed for
##               exact reproduction of your results.
##
## Return values:
## pcaRes    - pcaRes is the standart return class for all PCA methods
##          in this package.
##          Fields are:
##            pcaRes@completeObs - estimated complete observations
##            pcaRes@scores      - estimated scores
##            pcaRes@loadings    - estimated loadings (eigenvectors)
##            pcaRes@R2          - the individual R2 values
##            pcaRes@R2cum       - the cumulative R2 values
##            pcaRes@sDev        - standart deviation of the scores
##            pcaRes@nObs        - number of observations (rows of input matrix)
##            pcaRes@nVar        - number of variables (columns of input matrix)
##            pcaRes@centered    - boolean, TRUE if centered, FALSE otherwise
##            pcaRes@center      - row wise mean
##            pcaRes@subset      - subset used for estimation. Rows / columns
##                                 that contain only NAs are left out.
##            pcaRes@varLimit    - NULL, not used here
##            pcaRes@scaled      - boolean, TRUE if data was scaled, FALSE otherwise
##            pcaRes@nPcs        - number of principal components (d)
##            pcaRes@method      - "ppca"
##            pcaRes@missing     - Number of NAs in the data
##
## Requires:     MASS
##
## Author:    Wolfram Stacklies
##            Max Planck Institut fuer Molekulare Pflanzenphysiologie
##            Golm, Germany
## Date:      13.04.2006
##
## Contact:    wolfram.stacklies@gmail.com
##
##########################################################################################
          
ppca <- function(Matrix, nPcs = 2, center = TRUE, completeObs = TRUE, seed = NA, ... ) {

    ## Set the seed to the user defined value. This affects the generation
    ## of random values for the initial setup of the loading matrix
    if (!is.na(seed)) 
        set.seed(seed)

    d <- nPcs

    ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix)
    ## And now check if everything is right...
    if ( !checkData(Matrix, verbose = interactive()) ) {
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")
    }

    threshold <- 1e-5

    N <- nrow(Matrix)
    D <- ncol(Matrix)

    if (d > D) {
        stop("more components than matrix columns selected, exiting.\n")
    }

    Obs <- !is.na(Matrix)
    hidden <- which(is.na(Matrix))
    missing <- length(hidden)

    ## compute data mean and center data
    M <- NULL
    if(missing) {
        for (i in 1:D) {
            M[i] <- mean( Matrix[Obs[,i] ,i] )
        }
    } else { 
        M <- apply(Matrix, 2, mean) 
    }

    if (center) {
        Ye <- Matrix - repmat(M, N, 1)
    } else
        Ye <- Matrix

    if(missing) { Ye[hidden] <- 0 } 

## ------- Initialization
    r     <- sample(N)
    C    <- t(Ye[r[1:d], ,drop = FALSE])
    ## Random matrix with the same dimnames as Ye
    C    <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
    CtC    <- t(C) %*% C
    ## inv(C'C) C' X is the solution to the EM problem
    X    <- Ye %*% C %*% solve(CtC)
    recon    <- X %*% t(C)
    recon[hidden] <- 0
    ss    <- sum(sum((recon - Ye)^2)) / (N - missing)

    count <- 1
    old <- Inf
    
## ------ EM iterations
    while (count > 0) {
        ## E-step, (co)variances
        Sx <- solve(diag(d) + CtC/ss) 
        ss_old <- ss
        if(missing) {
            proj <- X %*% t(C)
            Ye[hidden] <- proj[hidden]
        }
        
        ## E step: expected values
        X <- Ye %*% C %*% Sx / ss
    
        ## M-step
        SumXtX <- t(X) %*% X

##        Replace the right matrix division from matlab
        C <- (t(Ye) %*% X) %*% solve( (SumXtX + N * Sx) )
        
        CtC <- t(C) %*% C
        ss <- ( sum(sum( (C %*% t(X) - t(Ye))^2 )) + N * sum(sum(CtC %*% Sx)) +
            missing * ss_old ) / (N * D)

        ## Some of the values may be negative at the beginning of the iteration,
        ## check that we are ot trying to calculate a log(<0)
        if( (ss < 0) | (det(Sx) < 0) | (ss_old < 0) ) {
            objective <- NaN
        } else {
            objective <- N * (D * log(ss) + sum(diag(Sx)) - log(det(Sx)) ) +
                sum(diag(SumXtX)) - missing * log(ss_old)
        }

        rel_ch <- abs( 1 - objective / old )
        old <- objective

        count <- count + 1
        if (!is.nan(rel_ch))  {
            if( (rel_ch < threshold) && (count > 5) ) {
                count <- 0
            }
        } else if (count > 1000) {
            count <- 0
            warning("Stopped after 1000 iterations, but rel_ch was NaN\n",
                "Results may be inaccurate\n")
        }    
    } ## End EM iteration
    C <- orth(C)
    evs <- eigen( cov(Ye %*% C) )
    vals <- evs[[1]]
    vecs <- evs[[2]]
    
    C <- C %*% vecs
    X <- Ye %*% C

    ## add data mean to expected complete data
    if (center)
        Ye <- Ye + repmat(M,N,1)

    ## Paramters in original Matlab implementation were:
    ## C (D by d)    - C has the approximate loadings (eigenvectors of the covariance matrix)
    ##          as columns.
    ## X        - The approximate scores 
    ## Ye (N by D)    - Expected complete observations.
    ## M (D by 1)    - Column wise data mean
    ## ss (scalar)    - isotropic variance outside subspace

    ## Replace NA values by estimated complete observations
    cObs <- Matrix
    cObs[hidden] <- Ye[hidden]

        ## Calculate R2cum
    R2cum <- NULL
    scaled  <- scale(cObs, center = TRUE, scale = FALSE)
        for (i in 1:ncol(C)) {
                difference <- scaled - ( X[,1:i, drop=FALSE] %*% t(C[,1:i, drop=FALSE]) )
                R2cum <- cbind( R2cum, 1 - ( sum(difference^2) / sum(scaled^2) ) )
        }

    ## Calculate R2
        R2 <- vector(length=length(R2cum), mode="numeric")
        R2[1] <- R2cum[1]
        if (ncol(C) > 1) {
                for (i in 2:ncol(C)) {
                        R2[i] <- R2cum[i] - R2cum[i - 1]
                }
        }


    ####################################################################
    ## Store the values in the pcaRes class, for compatibilty with
    ## other methods provided in the package
    ####################################################################

    result                 <- new("pcaRes")

    if (completeObs)
        result@completeObs <- cObs
    result@center          <- attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center")

    result@centered        <- center
    result@scaled          <- "none"
    result@scores          <- X
    result@loadings        <- C
    result@R2cum           <- c(R2cum)
    result@R2              <- R2
    result@sDev            <- apply(X, 2, sd)
    result@nObs            <- nrow(Matrix)
    result@nVar            <- ncol(Matrix)
##    result@subset        <- NULL
##    result@varLimit      <- NULL
    result@nPcs            <- ncol(C)
    result@method          <- "ppca"
    result@missing         <- sum(is.na(Matrix))
    
    return(result)
    
}

