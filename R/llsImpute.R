##' Missing value estimation using local least squares (LLS).  First,
##' k variables (for Microarrya data usually the genes)  are selected
##' by pearson, spearman or kendall correlation coefficients.  Then
##' missing values are imputed by a linear combination of the k
##' selected variables. The optimal combination is found by LLS
##' regression.  The method was first described by Kim et al,
##' Bioinformatics, 21(2),2005.
##'
##' Missing values are denoted as \code{NA}\cr It is not recommended
##' to use this function directely but rather to use the nni() wrapper
##' function. The methods provides two ways for missing value
##' estimation, selected by the \code{allVariables} option. The first
##' one is to use only complete variables for the  regression. This is
##' preferable when the number of incomplete variables is relatively
##' small.
##' 
##' The second way is to consider all variables as candidates for the
##' regression.  Hereby missing values are initially replaced by the
##' columns wise mean.  The method then iterates, using the current
##' estimate as input for the regression until the change between new
##' and old estimate falls below a threshold (0.001).
##' 
##' @title LLSimpute algorithm
##' @param Matrix \code{matrix} -- Data containing the variables
##' (genes) in columns and observations (samples) in rows. The data
##' may contain missing values, denoted as \code{NA}.
##' @param k \code{numeric} -- Cluster size, this is the number of
##' similar genes used for regression.
##' @param center \code{boolean} -- Mean center the data if TRUE
##' @param completeObs \code{boolean} -- Return the estimated complete
##' observations if  TRUE. This is the input data with NA values
##' replaced by the estimated values.
##' @param correlation \code{character} -- How to calculate the
##' distance between genes.  One out of pearson | kendall | spearman ,
##' see also help("cor").
##' @param allVariables \code{boolean} -- Use only complete genes to
##' do the regression if TRUE, all genes if FALSE.
##' @param maxSteps \code{numeric} -- Maximum number of iteration
##' steps if allGenes = TRUE.
##' @param xval \code{numeric} Use LLSimpute for cross
##' validation. xval is the index of the gene to estimate, all other
##' incomplete genes will be ignored if this parameter is set. We do
##' not consider them in the cross-validation.
##' @param verbose \code{boolean} -- Print step number and relative
##' change if TRUE and  allVariables = TRUE
##' @param ... Reserved for parameters used in future version of the
##' algorithm
##' @note Each step the generalized inverse of a \code{miss} x k
##' matrix is calculated. Where \code{miss} is the number of missing
##' values in  variable j and \code{k} the number of neighbours. This
##' may be slow for large values of k and / or many missing
##' values. See also help("ginv").
##' @return   \item{nniRes}{Standard nni (nearest neighbour
##' imputation) result object of this package. See
##' \code{\link{nniRes}} for details.}
##' @seealso \code{\link{pca}, \link{nniRes}, \link{nni}}.
##' @examples
##' ## Load a sample metabolite dataset (metaboliteData) with already 5\% of
##' ## data missing
##' data(metaboliteData)
##' ## Perform llsImpute using k = 10
##' ## Set allVariables TRUE because there are very few complete variables
##' result <- llsImpute(metaboliteData, k = 10, correlation="pearson", allVariables=TRUE)
##' ## Get the estimated complete observations
##' cObs <- completeObs(result)
##' @keywords multivariate
##' @export
##' @references Kim, H. and Golub, G.H. and Park, H.  - Missing value
##' estimation for DNA microarray gene expression data: local least
##' squares imputation.  \emph{Bioinformatics, 2005; 21(2):187-198.}
##' 
##' Troyanskaya O. and Cantor M. and Sherlock G. and Brown P. and
##' Hastie T. and Tibshirani R. and Botstein D. and Altman RB.  -
##' Missing value estimation methods for DNA microarrays.
##' \emph{Bioinformatics. 2001 Jun;17(6):520-525.}
##' @author Wolfram Stacklies
llsImpute <- function(Matrix, k=10, center=FALSE, completeObs=TRUE,
                      correlation="pearson", 
                      allVariables=FALSE, maxSteps=100, xval=NULL,
                      verbose=FALSE, ...) {

    threshold <- 0.001

    correlation <- match.arg(correlation, c("pearson", "kendall", "spearman"))

    ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix, rownames.force=TRUE)
    ## And now check if everything is right...
    if ( !checkData(Matrix, verbose = interactive()) ) {
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")
    }

    ## Exit if number of neighbours exceeds number of columns
    if (k > ncol(Matrix))
        stop("Cluster size larger than the number of columns, choose a k < ncol(Matrix)!")
 
    ## Set allVariables TRUE if k exceeds number of complete genes
    ## Print warning messages in the first case and when less than 50% of all genes are complete
    ## and allVariables == FALSE
    cg <- sum( apply(is.na(Matrix), 2, sum) == 0)
    if ( (k > cg) && (!allVariables) ) {
        warning("Cluster size larger than number of complete genes, using allVariables = TRUE")
        allVariables <- TRUE
    } else if ( (cg < (ncol(Matrix) / 2)) && (!allVariables) ) {
        warning("Less than 50% of the genes are complete, consider using allVariables = TRUE")
    } else if (sum(is.na(Matrix)) == 0)
        stop("No missing values, no need for missing value imputation :))")

    ## Find all genes with missing values
    missing <- apply(is.na(Matrix), 2, sum) > 0
    missIx <- which(missing == TRUE)
    # For cross validation we want to only estimate one variable, the others
    # are not considered in the cross validation anyway
    if (!is.null(xval))
        missIx = xval
    obs <- Matrix    ## working copy of the data
    Ye <- Matrix     ## Estimated complete observations

    ## Center the data column wise
    if (center) {
        obs   <- scale(Matrix, center = TRUE, scale = FALSE)
        Ye    <- obs
        means <- attr(Ye, "scaled:center")
    }

    if (allVariables) {
        compIx <- 1:ncol(obs)
        ## Impute the row average
        rowMeans <- apply(obs, 1, mean, na.rm = TRUE)
        for (i in 1:nrow(obs)) {
            obs[i, is.na(Matrix[i,])] <- rowMeans[i]
        }
        ## distances between all genes, ignore the diagonal (correlation to itself)
        distance = abs(cor(obs, obs, method = correlation))
    } else {
        compIx <- which(missing == FALSE)
        ## missing genes are the rows, complete genes the columns
        distance = abs(cor(obs[,missIx, drop=FALSE], obs[,compIx, drop=FALSE], use="pairwise.complete.obs",
                       method = correlation))
    }

    change <- Inf
    step <- 0
    while ( (change > threshold) && (step < maxSteps) ) {
        step <- step + 1
        iteration <- 0
        
        ## Do the regression and imputation
        for (index in missIx) {
            iteration <- iteration + 1
	    if (allVariables) {
                similar <- sort(distance[iteration,], index.return = TRUE, decreasing = TRUE)
                simIx <- compIx[ similar$ix[similar$ix != iteration][1:k] ]
            } else {
                similar <- sort(distance[iteration,], index.return = TRUE, decreasing = TRUE)
                simIx <- compIx[ similar$ix[1:k] ]
            }

            ##
            ## Do a regression against the k most similar genes
            ## See Kim et. al 2005 for details
            ##
            target <- obs[, index, drop = FALSE]
            tMiss <- is.na(Matrix[, index, drop = FALSE])

            Apart <- obs[!tMiss, simIx, drop = FALSE]
            Bpart <- obs[tMiss, simIx, drop = FALSE]
            targetComplete <- target[!tMiss, , drop = FALSE]
            X <- MASS::ginv(Apart) %*% targetComplete
            estimate <- Bpart %*% X

            ## Impute the estimate
            Ye[tMiss, index] <- estimate
        }

        ## We do not want to iterate if allVariables == FALSE
        if (!allVariables || !is.null(xval)) {
            break
        } else {
            ## relative change in estimation
            change <- sqrt(sum( (obs - Ye)^2 ) / sum(obs^2))
            obs <- Ye
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
    result <- new("nniRes")

    if(completeObs) {
        Ye[!is.na(Matrix)] <- Matrix[!is.na(Matrix)]
        result@completeObs <- Ye
    }
    result@centered        <- center
    result@center          <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
    result@nObs            <- nrow(Matrix)
    result@nVar            <- ncol(Matrix)
    result@method          <- "llsImpute"
    result@correlation     <- correlation
    result@k               <- k
    result@missing         <- sum(is.na(Matrix))

    return(result)        
} 

