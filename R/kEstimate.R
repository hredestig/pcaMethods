##########################################################################################
##
## Perform cross validation to estimate an optimal number of components for missing
## value estimation.
## Cross validation is done on the complete subset of a varialbe (gene).
## For each incomplete gene, the available values are diveded into a user defined
## number of segments. The segments have equal size, but are chosen from a random
## equal distribution. The available values are covered completely.
## PPCA, BPCA, SVDimpute, Nipals PCA and llsImpute may be used for imputation.
##
## The whole cross validation is repeated several times.
## As error measure the NRMSEP (see Feten 2005) is used. This error basically
## normalizes the RMSD between original data and estimate by the gene-wise
## varicance. This makes sense because a higher variance will lead to a
## higher estimation error.
##
## Parameters:
##    Matrix      - numeric matrix containing observations in rows and 
##                  genes in columns
##    method      - One of ppca | bpca | svdImpute | nipals | llsImpute | llsImputeAll | nlpca.
##                  llsImputeAll uses llsImpute with option allVariables = TRUE.
##    evalPcs     - The principal components to use for cross validation
##                  or cluster sizes if used with llsImpute.
##                  Should be an array containing integer values, eg. evalPcs = 1:10
##                  or evalPcs = C(2,5,8).
##                  The NRMSEP is calculated for each component.
##    segs        - number of segments for cross validation
##    nruncv      - Times the whole cross validation is repeated
##    em          - Error measure (em) to use, this can be "nrmsep" (default)
##                  or "q2"
##    allVariables- If TRUE, the NRMSEP is calculated for all genes,
##                  If FALSE, only the incomplete ones are included.
##                  You maybe want to do this to compare several methods on a 
##                  complete data set
##    verbose     - If TRUE, the NRMSEP and the variance are prited. 
##    random      - Impute normal distributed random values with
##                  same mean and standard deviation that the origial data.
##                  This is only thought for comparison.
##
## Return values:
## A list with the elements
## mink      - number of PCs for which the minimal average NRMSEP was obtained
## eError    - an array of of size length(evalPcs). Contains the average error
##             of the cross validation runs for each number of components.
## eErrorMat - Estimation error for each of the separately done cross validations.
##             Rows represent the separate runs, columns the number of components
##             used for estimation.
## evalPcs   - The evaluated numbers of pc's or cluster sizes (same as input)
##
## Author:  Wolfram Stacklies
##          Max Planck Institute for Molecular Plant Physiology
## Date:    06/28/2006
## Contact: wolfram.stacklies@gmail.com
##
## Last modified: 05/02/07 by Wolfram Stacklies
##
##########################################################################################

kEstimate <- function(Matrix, method = "ppca", evalPcs = 1:3, segs = 3, nruncv = 10,
                      em = "nrmsep", allVariables = FALSE, verbose = interactive(),
                      random = FALSE,...) {

    method <- match.arg(method, c("ppca", "bpca", "svdImpute", "nipals", 
                                  "llsImpute", "llsImputeAll", "nlpca"))
    em <- match.arg(em, c("nrmsep", "q2"))
    maxPcs <- max(evalPcs)
    lengthPcs <- length(evalPcs)

     ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix)
    if(maxPcs > (ncol(Matrix) - 1))
        stop("maxPcs exceeds matrix size, choose a lower value!")

     ## And now check if everything is right...
    if( !checkData(Matrix, verbose=interactive()) )
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

    if( (sum(is.na(Matrix)) == 0) && (allVariables == FALSE) )
        stop("No missing values. Maybe you want to set allVariables = TRUE. Exiting\n")


    missing <- apply(is.na(Matrix), 2, sum) > 0
    missIx     <- which(missing == TRUE)
    if (allVariables)
        missIx <- 1:ncol(Matrix)

    complete <- !missing
    compIx    <- which(complete == TRUE)

    finalError <- matrix(0, nruncv, lengthPcs)
    iteration <- 0
    for(nPcs in evalPcs) {
        iteration = iteration + 1
        if (verbose) { cat("Doing CV for ", nPcs, " component(s) ") }
        error <- matrix(0, length(missIx), nruncv)
        for(cviter in 1:nruncv) {
            if (verbose) { cat(".") }
            pos <- 0
            for (index in missIx) {
                pos <- pos + 1
                target <- Matrix[, index, drop = FALSE]
                compObs <- !is.na(target)
                missObs <- is.na(target)
                nObs <- sum(compObs)

                ## Remove all observations that are missing in the target genes,
                ## as additional missing values may tamper the results
                set <- Matrix[compObs,]

                if (nObs >= (2 * segs)) {
                    segments <- segs
                } else
                    segments <- ceiling(nObs / 2)

                ## We assume normal distributed missing values when choosing the segments
                cvsegs <- cvsegments(nObs, segments)
                set <- Matrix[compObs,]
                nrmsep <- 0
                q2 <- 0
    
                for (i in length(cvsegs)) {
                    n <- length(cvsegs[[i]]) ## n is the number of created missing values

                    ## Impute equally distributed random numbers, for testing only
                    if (random) {
                        original <- set[cvsegs[[i]], index]
                        estimate <- rnorm(n, mean(original), sd(original))
                    } else {
                    ## Impute values using the given regression method
                        testSet <- set
                        testSet[cvsegs[[i]], index] <- NA
                        if (method == "llsImpute") {
                            estimate <- llsImpute(testSet, k = nPcs, verbose = FALSE,
                                                  allVariables = FALSE, 
                                                  center = FALSE)@completeObs
                        } else if (method == "llsImputeAll") {
                            estimate <- llsImpute(testSet, k = nPcs, verbose = FALSE,
                                                  allVariables = TRUE, 
                                                  center = FALSE)@completeObs
                        } else {
                            estimate <- pca(testSet, nPcs = nPcs, verbose = FALSE,
                                            method = method, center = TRUE,...)@completeObs
                        }
                        estimate <- estimate[, index]
                        original <- target[compObs, ]
                    }
                    ## Error of prediction, error is calculated for removed elements only
                    if (em == "nrmsep") {
                        nrmsep <- nrmsep +
                            sum( (original - estimate)^2) /
                            (n * sum( (original - mean(original))^2)  / (nObs - 1) )
                    } else {
                        q2 <- q2 +
                              sum( (original - estimate)^2 ) / sum(original^2)
                    }
                } ## iteration over cv segments
        
                if (em == "nrmsep") {
                    error[pos, cviter] <- (nrmsep / length(cvsegs))
                } else
                    error[pos, cviter] <- 1 - (q2 / length(cvsegs))
            } ## iteration over genes
        } ##iteration over nruncv
        
        errorAllGenes <- apply(error, 2, sum) / nrow(error)
        if (em == "nrmsep")
            finalError[, iteration] <- sqrt(errorAllGenes)
        else
            finalError[, iteration] <- errorAllGenes

        if (verbose)
            cat(" - the average", em, "is", 
                sum(finalError[,iteration]) / nrow(finalError), 
                ". The variance is", var(finalError[,iteration]), "\n")
    } ## iteration over number components

    avgError <- apply(finalError, 2, sum) / nrow(finalError)
    ret <- list()
    ret$mink <- evalPcs[which(avgError == min(avgError))]
    ret$eError <- avgError
    ret$eErrorMat <- finalError
    ret$evalPcs <- evalPcs

    return(ret)
}

