#####################################################################################
##
## svdImpute <- function(Matrix, nPcs = 5, center = TRUE, completeObs = TRUE, threshold = 0.01, 
##                       maxSteps = 100, verbose = interactive(), ...)
##
## Implements the SVDimpute algorithm as described in Troyanskaya 2001.
## Initially all missing values are replaced with 0. The
## algorithm then selects the 'nPcs' most significant Eigengenes.
## Eigengenes are the eigenvectors when genes are considered as samples.
## Missing values are estimated by regressing the target gene against
## the 'pPcs' most significant Eigengenes.
##
## Parameters:
## Matrix      - A numeric data matrix. Missing values are 
##               denoted as 'NA'
## nPcs        - Number of components used for estimation
## center      - mean center the data if TRUE
## completeObs - Return estimated complete observations if TRUE
## threshold   - The iteration stops if the change in the matrix
##               falls below this threshold. Default is 0.01.
##               (0.01 was empirically determined by 
##               Troyanskaya et. al)
## maxSteps    - Maximum number of iteration steps. Default is 100
## verbose     - Print some output if TRUE. Default is interactive()
##
## Author:       Wolfram Stacklies
##               Max Planck Institut fuer Molekulare Pflanzenphysiologie
##               Golm, Germany
## Date:         08/06/2006
##
## Contact:      wolfram.stacklies@gmail.com
##
#####################################################################################

svdImpute <- function(Matrix, nPcs = 2, center = TRUE, completeObs = TRUE, threshold = 0.01, 
                      maxSteps = 100, verbose = interactive(), ...) {

    Matrix <- as.matrix(Matrix)

    if (!checkData(Matrix, verbose = verbose))
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

    if (nPcs > ncol(Matrix))
        stop("more components than matrix columns selected, exiting")

    Ye <- Matrix
    missing <- is.na(Matrix)
    temp <- apply(missing, 2, sum)
    missIx <- which(temp != 0)

    ## Center the data column wise
    if (center) {
        Ye <- scale(Matrix, center = TRUE, scale = FALSE)
	  means <- attr(Ye, "scaled:center")
    }

    ## Initially set estimates to 0
    Ye[missing] <- 0

    ## Now do the regression
    count <- 0
    error <- Inf

    while ( (error > threshold) && (count < maxSteps) ) {
        res         <- prcomp(t(Ye), center = FALSE, scale = FALSE, retx = TRUE)
        loadings    <- res$rotation[,1:nPcs, drop = FALSE]
        sDev        <- res$sdev

        ## Estimate missing values as a linear combination of the eigenvectors
        ## The optimal solution is found by regression against the k eigengenes
        for (index in missIx) {
            target <- Ye[!missing[,index],index, drop = FALSE]
            Apart <- loadings[!missing[,index], , drop = FALSE]
            Bpart <- loadings[missing[,index], , drop = FALSE]
            X <- ginv(Apart) %*% target
            estimate <- Bpart %*% X
        
            Ye[missing[,index], index] <- estimate
        }
            
        count <- count + 1
        if (count > 5) {
            error <- sqrt(sum( (YeOld - Ye)^2 ) / sum(YeOld^2))
            if (verbose) { cat("change in estimate: ", error, "\n") }
        }
        YeOld <- Ye
    }

    ## Add the original mean
    if (center) {
        for(i in 1:ncol(Ye)) {
            Ye[,i] <- Ye[,i] + means[i]
        }
    }

    ## Calculate R2cum
    tmp <- prcomp(Ye, center = FALSE, scale = FALSE, retx = TRUE)
    loadings <- tmp$rotation
    scores <- tmp$x

    R2cum <- NULL
    centered  <- scale(Ye, center = TRUE, scale = FALSE)
    for (i in 1:nPcs) {
        difference <- centered - ( scores[,1:i, drop=FALSE] %*% t(loadings[,1:i, drop=FALSE]) )
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


    ########################################################################
    ## Store the values in the pcaRes class, for compatibilty with
    ## other methods provided in the package
    ########################################################################
    
    result <- new("pcaRes")

    cObs <- Matrix
    cObs[missing] <- Ye[missing]
    if (completeObs)
        result@completeObs <- cObs
    result@center <- attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center") 
    result@centered <- center
    result@scaled <- "none"
    result@scores <- scores
    result@loadings <- loadings
    result@R2cum <- c(R2cum)
    result@R2 <- R2
    result@sDev <- sDev[1:nPcs]
    result@nObs <- nrow(Matrix)
    result@nVar <- ncol(Matrix)
    result@nPcs <- nPcs
    result@method <- "svdImpute"
    result@missing <- sum(is.na(Matrix))

    return(result)

}
