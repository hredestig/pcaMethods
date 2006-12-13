####################################################################################################
## LLSimpute <- function(Matrix, k = 10, center = TRUE, completeObs = TRUE, correlation = "pearson", 
##                      allGenes = FALSE, maxSteps = 100, verbose = interactive(), ...) {
##
## Missing value estimation using local least sqares (LLS).
## First k genes are selected by pearson or spearman correlation coefficients.
## Then missing values are imputed by regression against the k selected
## genes.
## The method was first described in Kim et al, Bioinformatics, 21(2),2005.
## The allGenes option allows to choose if either only complete genes or all
## genes may be used for estimation. If all genes are used, initially the missing
## values are replaced by the columns wise mean. The method then iterates, using
## the current estimate as input for the regression until the change between new
## and old estimate falls below a threshold.
##
## The method considers columns as variables (genes) and rows as 
## observations (samples)
##
## Parameters:
## Matrix      - A numeric matrix or data frame. Missing values are denoted
##               as 'NA'
## k           - Cluster size, the number of similar genes used for estimation
## center      - mean center the data if TRUE
## completeObs - Return estimated complete observations if TRUE
## correlation - one out of "pearson | kendall | spearman". See also help("cor").
##               allGenes If FALSE only complete genes are used for the regression, if set
##               TRUE all genes are considered. Therefore missing values are initially replaced 
##               by the row wise mean. Then the estimation is repeated with the newly imputed
##               values until the change falls below a certain threshold (here 0.001).
## maxSteps    - Maximum number of steps when allGenes = TRUE
## verbose     - Print step number and relative change if TRUE and allGenes = TRUE
##
## Return values:
## nniRes      - a nearest neighbour imputation (nni) result object
##
## Author:    Wolfram Stacklies
##            MPG/CAS Partner Institute for Computational Biology (PICB)
##            Shanghai, P.R. China
## Date:      12/11/2006
##
## Contact:    wolfram.stacklies@gmail.com
##
####################################################################################################

llsImpute <- function(Matrix, k = 10, center = FALSE, completeObs = TRUE, correlation = "pearson", 
                      allGenes = FALSE, maxSteps = 100, verbose = interactive(), ...) {

    threshold = 0.001

    correlation = match.arg(correlation, c("pearson", "kendall", "spearman"))

    ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix)
    ## And now check if everything is right...
    if ( !checkData(Matrix, verbose = interactive()) ) {
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")
    }

    ## Exit if number of neighbours exceeds number of columns
    if (k > ncol(Matrix))
        stop("Cluster size larger than the number of columns, choose a k < ncol(Matrix)!")
 
    ## Set allGenes TRUE if k exceeds number of complete genes
    ## Print warning messages in the first case and when less than 50% of all genes are complete
    ## and allGenes == FALSE
    cg = sum( apply(is.na(Matrix), 2, sum) == 0)
    if ( (k > cg) && (!allGenes) ) {
        warning("Cluster size larger than number of complete genes, using allGenes = TRUE")
        allGenes = TRUE
    } else if ( (cg < (ncol(Matrix) / 2)) && (!allGenes) ) {
        warning("Less than 50% of the genes are complete, consider using allGenes = TRUE")
    } else if (sum(is.na(Matrix)) == 0)
        stop("No missing values, no need for missing value imputation :))")

    ## Find all genes with missing values
    missing = apply(is.na(Matrix), 2, sum) > 0
    missIx = which(missing == TRUE)
    obs = Matrix    ## working copy of the data
    Ye = Matrix     ## Estimated complete observations

    ## Center the data column wise
    if (center) {
        obs   <- scale(Matrix, center = TRUE, scale = FALSE)
        Ye    <- obs
        means <- attr(Ye, "scaled:center")
    }

    if (allGenes) {
        compIx = 1:ncol(obs)
        ## Impute the row average
        rowMeans = apply(obs, 1, mean, na.rm = TRUE)
        for (i in 1:nrow(obs)) {
            obs[i, is.na(Matrix[i,])] = rowMeans[i]
        }
        ## distances between all genes, ignore the diagonal (correlation to itself)
        distance = abs(cor(obs, obs, method = correlation))
    } else {
        compIx = which(missing == FALSE)
        ## missing genes are the rows, complete genes the columns
        distance = abs(cor(obs[,missIx], obs[,compIx], use="pairwise.complete.obs",
                       method = correlation))
    }

    change = Inf
    step = 0
    while ( (change > threshold) && (step < maxSteps) ) {
        step = step + 1
        iteration = 0
        
        ## Do the regression and imputation
        for (index in missIx) {
            iteration = iteration + 1
	    if (allGenes) {
                similar = sort(distance[iteration,], index.return = TRUE, decreasing = TRUE)
                simIx = compIx[ similar$ix[similar$ix != iteration][1:k] ]
            } else {
                similar = sort(distance[iteration,], index.return = TRUE, decreasing = TRUE)
                simIx = compIx[ similar$ix[1:k] ]
            }

            ##
            ## Do a regression against the k most similar genes
            ## See Kim et. al 2005 for details
            ##
            target = obs[, index, drop = FALSE]
            tMiss = is.na(Matrix[, index, drop = FALSE])

            Apart = obs[!tMiss, simIx, drop = FALSE]
            Bpart = obs[tMiss, simIx, drop = FALSE]
            targetComplete = target[!tMiss, , drop = FALSE]
            X = ginv(Apart) %*% targetComplete
            estimate = Bpart %*% X

            ## Impute the estimate
            Ye[tMiss, index] = estimate
        }

        ## We do not want to iterate if allGenes = FALSE
        if (!allGenes) {
            break
        } else {
            ## relative change in estimation
            change = sqrt(sum( (obs - Ye)^2 ) / sum(obs^2))
            obs = Ye
            if (verbose) {
                cat("Step number     : ", step, '\n')
                cat("Relative change : ", change, '\n')
                cat("---------------", '\n')
            }
        }
    }

    ## Add the original mean
    if (center) {
        for(i in 1:ncol(Ye)) {
            Ye[,i] <- Ye[,i] + means[i]
        }
    }

    ## Build the nniRes object
    ##
    result = new("nniRes")

    if(completeObs) {
        Ye[!is.na(Matrix)] = Matrix[!is.na(Matrix)]
        result@completeObs = Ye
    }
    result@center          = attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center")
    result@centered        = center
    result@scaled          = "none"
    result@nObs            = nrow(Matrix)
    result@nVar            = ncol(Matrix)
    result@method          = "llsImpute"
    result@correlation     = "correlation"
    result@k               = k
    result@missing         = sum(is.na(Matrix))

    return(result)        
} 

