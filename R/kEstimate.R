##########################################################################################
##
## Perform cross validation to estimate an optimal number of components for missing
## value estimation.
## Cross validation is done on the complete subset of a varialbe (gene).
## For each incomplete gene, the available values are diveded into a user defined
## number of segments. The segments have equal size, but are chosen from a random
## equal distribution. The available values are covered completely.
## PPCA, BPCA, SVDimpute, NLPCA, Nipals PCA and llsImpute may be used for imputation.
##
## The whole cross validation is repeated several times.
## As error measure the NRMSEP (see Feten 2005) or Q2 is used. The NRMSEP basically
## normalizes the RMSD between original data and estimate by the gene-wise
## varicance. This makes sense because a higher variance will lead to a
## higher estimation error.
##
## Parameters:
##    Matrix      - numeric matrix containing observations in rows and 
##                  genes in columns
##    method      - One of ppca | bpca | svdImpute | nipals | llsImpute | llsImputeAll | nlpca.
##                  llsImputeAll calls llsImpute with the allVariables = TRUE option.
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
##
## Return values:
## A list with the elements
## bestNPcs          - number of PCs for which the minimal average NRMSEP or the maximal 
##                     average Q2 was found. For llsImpute this returns the optimal k.
## eError            - an array of of size length(evalPcs). Contains the average error
##                     of the cross validation runs for each number of components.
## variableWiseError - Matrix of size(incomplete_variables x length(evalPcs)).
##                     This contains the NRMSEP or Q2 distance for each variable and each number of PC's.
##                     This allows to easily see for wich variables imputation makes sense and for
##                     which one it should not be done or mean imputation should be used.
## evalPcs           - The evaluated numbers of pc's or number of neighbours (same as input)
## variableIx        - Index of the incomplete variables. This can be used to map the variable
##                     wise error to the original data.
##
## Author:  Wolfram Stacklies
##          MPG-CAS Partner Institute for Computational Biology
## Date:    06/28/2006
## Contact: wolfram.stacklies@gmail.com
##
## Last modified: 14/02/07 by Wolfram Stacklies
##
##########################################################################################

kEstimate <- function(Matrix, method = "ppca", evalPcs = 1:3, segs = 3, nruncv = 5,
                      em = "q2", allVariables = FALSE, verbose = interactive(), ...) {

   fastKE <- FALSE
   if (method == "ppca" || method == "bpca" || method == "nipals" || method == "nlpca")
       fastKE <- TRUE

    method <- match.arg(method, c("ppca", "bpca", "svdImpute", "nipals", 
                                  "llsImpute", "llsImputeAll", "nlpca"))
    em <- match.arg(em, c("nrmsep", "q2"))
    maxPcs <- max(evalPcs)
    lengthPcs <- length(evalPcs)

     ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix, rownames.force=TRUE)
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

    error <- matrix(0, length(missIx), length(evalPcs))
    iteration <- 0
    for(nPcs in evalPcs) {
        # If the estimated observations are just scores %*% t(loadings) we can calculate
        # all we need at once, this saves many iterations...
        if (fastKE) nPcs = maxPcs

        iteration = iteration + 1
        if (verbose && !fastKE) { cat("Doing CV for ", nPcs, " component(s) \n") }
        else if (verbose && fastKE) {cat("Doing CV ... \n")}
        for(cviter in 1:nruncv) {
            pos <- 0
            if (verbose) cat("Incomplete variable index: ")
            for (index in missIx) {
                pos <- pos + 1
                cat(pos, ":", sep="")
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
                if (fastKE) {
                    nrmsep <- array(0, length(evalPcs))
                    q2 <- array(0, length(evalPcs))
                } else {
                    nrmsep <- 0; q2 <- 0
                }
    
                for (i in 1:length(cvsegs)) {
                    n <- length(cvsegs[[i]]) ## n is the number of created missing values
                    ## Impute values using the given regression method
                        testSet <- set
                        testSet[cvsegs[[i]], index] <- NA
                        if (method == "llsImpute") {
                            estimate <- llsImpute(testSet, k = nPcs, verbose = FALSE,
                                                  allVariables = FALSE, 
                                                  center = FALSE, xval = index)
                        } else if (method == "llsImputeAll") {
                            estimate <- llsImpute(testSet, k = nPcs, verbose = FALSE,
                                                  allVariables = TRUE, 
                                                  center = FALSE, xval = index)
                        } else {
                            estimate <- pca(testSet, nPcs = nPcs, verbose = FALSE,
                                            method = method, center = TRUE,...)
                        }

                        if (fastKE) {
                            for (np in evalPcs) {
                                estiFitted <- fitted(estimate, data = NULL, nPcs = np)
                                estimateVec <- estiFitted[, index]
                                original <- target[compObs, ]
                                estimateVec[-cvsegs[[i]]] <- testSet[-cvsegs[[i]], index]
                                ## Error of prediction, error is calculated for removed elements only
                                nIx <- which(evalPcs == np) 
                                if (em == "nrmsep") {
                                    nrmsep[nIx] <- nrmsep[nIx] + sum( (original - estimateVec)^2) 
                                 } else {
                                     q2[nIx] <- q2[nIx] + sum( (original - estimateVec)^2 )
                                 }
                             }    
                         } else {
                             estimate <- estimate@completeObs[, index]
                             original <- target[compObs, ]
                             ## Error of prediction, error is calculated for removed elements only
                             if (em == "nrmsep") {
                                 nrmsep <- nrmsep + sum( (original - estimate)^2)
                             } else {
                                 q2 <- q2 + sum( (original - estimate)^2 )
                             }
                         }
                } ## iteration over cv segments
                
                if (fastKE) {
                    if (em == "nrmsep") {
                        error[pos, ] <- error[pos, ] + nrmsep / (nrow(set) * var(set[,index]))
                    } else
                        error[pos, ] <- error[pos, ] + (1 - (q2 / sum(set[, index]^2)))
                } else {
                    if (em == "nrmsep") {
                        error[pos, iteration] <- error[pos, iteration] + 
                                                 nrmsep / (nrow(set) * var(set[,index]))
                    } else
                        error[pos, iteration] <- error[pos, iteration] + (1 - (q2 / sum(set[, index]^2)))
                }
            } ## iteration over variables
            if (verbose) cat("\n")
            
        } ##iteration over nruncv
        
        # The error is the sum over the independent cross validation runs
        error <- error / nruncv
        
        if (verbose && !fastKE)
            cat("The average", em, "for k =", iteration, "is", 
                sum(error[,iteration]) / nrow(error), "\n")

       # if nlpca, ppca, bpca, nipals we do not need to iterate over the number of components...
       if (fastKE) break
    } ## iteration over number components

    if (em == "nrmsep")
        avgError <- sqrt(apply(error, 2, sum) / nrow(error))
    else
        avgError <- apply(error, 2, sum) / nrow(error)

    ret <- list()
    if (em == "nrmsep") ret$bestNPcs <- evalPcs[which(avgError == min(avgError))]
    else ret$bestNPcs <- evalPcs[which(avgError == max(avgError))]
    ret$eError <- avgError
    if(em == "nrmsep") ret$variableWiseError <- sqrt(error)
    else ret$variableWiseError <- error
    ret$evalPcs <- evalPcs
    ret$variableIx <- missIx

    return(ret)
}
