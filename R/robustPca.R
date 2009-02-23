
#####################################################################################
##
## PCA using robustSVD.
## This is a PCA implementation robust to outliers in a data set. It
## can also handle missing values, it is however NOT intended to be
## used for missing value estimation.  As it is based on robustSVD we
## will get an accurate estimation for the loadings also on incomplete
## data or data with outliers.  If the data show missing values,
## scores are caluclated by just setting all NA - values to zero. This
## is not expected to produce accurate results. As scores are just
## calculated using Matrix %*% loadings they are of course affected by
## outliers.  Use one of the other methods coming with this package
## (like PPCA or BPCA) if you want to do missing value estimation.
##
## Parameters:
## Matrix      - A numeric matrix or data frame. Missing values are denoted
##               as 'NA'
## nPcs        - Number of principal components to calculate
## center      - Mean center the data, if TRUE
## completeObs - Return estimated complete observations if TRUE
##
## Return values:
## pcaRes - a pcaRes object. pcaRes is the standard return object for all PCA methods.
##
## Requires:  aroma.light
##
## Author:  Wolfram Stacklies
##          CAS-MPG Partner Institute for Computational Biology
## Date:    04/03/2007
## Contact: wolfram.stacklies@gmail.com
##
#####################################################################################

robustPca <- function(Matrix, nPcs = 2, center = TRUE, completeObs = FALSE, 
                      verbose = interactive(), ... ) {

    ## Do some basic checks
    Matrix <- as.matrix(Matrix, rownames.force=TRUE)
    if (!checkData(Matrix, verbose = verbose))
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

    if (nPcs > ncol(Matrix)) 
        stop("more components than matrix columns selected, exiting.\n")

    if (center) {
        mat <- scale(Matrix, center = TRUE, scale = FALSE)
        means <- attr(mat, "scaled:center")
    } else
        mat <- Matrix
    
    nas <- is.na(Matrix)
    complete <- sum(nas) == 0

    if (!complete && verbose) {
       cat("Input data is not complete.\n")
       cat("Scores, R2 and R2cum may be inaccurate, handle with care !!\n")
    }

    if (!complete && completeObs)
        warning("It is not recommended to use robustPca for missing value estimation!!\n")
    svdSol <- robustSvd(mat)
	
    ## Sort the eigenvalues and eigenvectors
    loadings <- svdSol$v[, 1:nPcs, drop = FALSE]
    sDev     <- svdSol$d[1:nPcs] / sqrt(max(1, nrow(Matrix) - 1))

    ## We estimate the scores by just setting all NA values to 0
    ## This is a bad approximation, I know... Use ppca / bpca or other missing value
    ## estimation methods included in this package
    compMat <- mat
    compMat[is.na(compMat)] <- 0
    scores   <- compMat %*% loadings

    ## Calculate R2cum (on the complete observations only)
    R2cum <- NULL
    centered  <- scale(Matrix, center = TRUE, scale = FALSE)
    for (i in 1:nPcs) {
        difference <- centered[!nas] - ( scores[,1:i, drop=FALSE] %*% t(loadings[,1:i, drop=FALSE]) )[!nas]
        R2cum <- cbind( R2cum, 1 - ( sum(difference^2) / sum(centered[!nas]^2) ) )
    }

    ## Calculate R2
    R2 <- vector(length=length(R2cum), mode="numeric")
    R2[1] <- R2cum[1]
    if (nPcs > 1) {
        for (i in 2:nPcs) {
            R2[i] <- R2cum[i] - R2cum[i - 1]
        }
    }

    result <- new("pcaRes")
    result@centered <- center
    result@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
    result@loadings <- loadings
        colnames(result@loadings) <- paste("PC", 1:nPcs, sep = "")
        rownames(result@loadings) <- colnames(Matrix)
    result@sDev <- sDev
    if (completeObs) {
        Ye <- scores %*% t(loadings)
        if (center) {
            for(i in 1:ncol(Ye))
                Ye[,i] <- Ye[,i] + means[i]
        }
        cObs <- Matrix
        cObs[is.na(Matrix)] <- Ye[is.na(Matrix)]
        result@completeObs <- cObs
    }
    result@scores <- scores
        colnames(result@scores) <- paste("PC", 1:nPcs, sep = "")
        rownames(result@scores) <- rownames(Matrix) 
    result@R2cum <- c(R2cum)
    result@R2 <- R2
    result@nObs <- nrow(Matrix)
    result@nVar <- ncol(Matrix)
    result@nPcs <- nPcs
    result@method <- "robustPca"
    result@missing <- sum(is.na(Matrix))

    return(result)
}


###############################################################################
##
## Singular Value Decomposition using alternating L1 norm.
## The method is robust to outliers.
##  
## See: Hawkins, Douglas M, Li Liu, and S Stanley Young (2001)
## Robust Singular Value Decomposition,
## National Institute of Statistical Sciences, Technical Report Number 122.
## http://www.niss.org/technicalreports/tr122.pdf
##
## This function needs the 'weightedMedian' function provided by the
## Bioconductor 'aroma.light' package. Visit www.biconductor.org for details.
##  
## If rank(x)<ncol(x) (nrow(x)<ncol(x)) this may fail...???
##
## Thanks a lot, Kevin!!!
##
## Author:	Kevin Wright.
## 		Adaptions for the pcaMethods package made by Wolfram Stacklies
## Date:	04/04/2007
###############################################################################

robustSvd <- function(x) {
	
    ## We need the weightedMedian function provided by the aroma.light
    ## package. However we do not want to make the whole package dependant
    ## on aroma.light
    if (!require(aroma.light, quietly = TRUE))
        stop("The aroma.light package is required in order to use this function.
The package is available at www.bioconductor.org")

    ## Do some basic checks
    x <- as.matrix(x)
    if (!checkData(x))
      stop("Invalid data format! Use checkData(x, verbose = TRUE) for details.\n")

    ## Define a couple of helper functions
    L1RegCoef <- function(x,a){
      keep <- (a!=0) & (!is.na(x))
      a <- a[keep]
    return ( weightedMedian(x[keep]/a,abs(a), interpolate = FALSE) )
    }

    L1Eigen <- function(x,a,b){
        x <- as.vector(x) # Convert from matrix to vector
        ab <- as.vector(outer(a,b))
        keep <- (ab!=0) & (!is.na(x))
        ab <- ab[keep]
        return( weightedMedian(x[keep]/ab,abs(ab), interpolate = FALSE) )
    }

    ## Initialize outputs
    svdu <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
    svdv <- matrix(NA,nrow=ncol(x),ncol=ncol(x))
    svdd <- rep(NA,ncol(x))

    for(k in 1:ncol(x)) {
        ak <- apply(abs(x),1,median,na.rm=TRUE)
        converged <- FALSE

        while(!converged) {
            akprev <- ak
            c <- apply(x,2,L1RegCoef,ak)
            bk <- c/sqrt(sum(c^2))
            d <- apply(x,1,L1RegCoef,bk)
            ak <- d/sqrt(sum(d^2))
            if(sum((ak-akprev)^2)< 1e-10) converged <- TRUE
        }
        eigenk <- L1Eigen(x,ak,bk)
        ## Deflate the x matrix
        x <- x - eigenk * ak %*% t(bk)
        ## Store eigen triple for output
        svdu[,k] <- ak
        svdv[,k] <- bk 
        svdd[k] <- eigenk
    }

    ## Create the result object
    ret <- list()

    ret$d <- svdd
    ret$u <- svdu
    ret$v <- svdv
    return(ret)
}

