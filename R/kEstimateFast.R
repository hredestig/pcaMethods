##########################################################################################
## 
## This is a simple estimator for the optimal number of componets / similar genes
## when applying pca or llsImpute for missing value estimation.
## No cross validation is performed, instead the estimation quality is defined
## as Matrix[!missing] - Estimate[!missing]. This will give a relatively rough
## estimate, but the number of iterations equals the length of parameter evalPcs.
## Does not work with LLSimpute.
##
## Parameters:
##    Matrix      - numeric matrix containing observations in rows and 
##                  genes in columns
##    method      - One of ppca | bpca | svdImpute | nipals | nlpca.
##    evalPcs     - The principal components to use for cross validation
##                  or cluster sizes if used with llsImpute.
##                  Should be an array containing integer values, eg. evalPcs = 1:10
##                  or evalPcs = C(2,5,8).
##                  The NRMSEP is calculated for each component.
##    em          - Error measure (em) to use, this can be "nrmsep" (default)
##                  or "q2"
##    verbose     - If TRUE, the NRMSEP and the variance are printed out. 
##
## Return values:
## A list with the elements
## mink      - number of PCs for which the minimal average NRMSEP was obtained
## eError    - an array of of size length(evalPcs). Contains the error found for 
##             each number of components.
## evalPcs   - The evaluated numbers of pc's or cluster sizes (same as input)
##
## Author:  Wolfram Stacklies
##          CAS-MPG Partner Institute for Computational Biology
## Date:    02/09/2007
## Contact: wolfram.stacklies@gmail.com
##
##########################################################################################

kEstimateFast <- function(Matrix, method = "ppca", evalPcs = 1:3,
                          em = "nrmsep", verbose = interactive(), ...) {

    method <- match.arg(method, c("ppca", "bpca", "svdImpute", "nipals", "nlpca"))
    em <- match.arg(em, c("nrmsep", "q2"))
    maxPcs <- max(evalPcs)
    lengthPcs <- length(evalPcs)
    missing <- is.na(Matrix)
    error <- array(0, lengthPcs)

     ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix)
    if(maxPcs > (ncol(Matrix) - 1))
        stop("maxPcs exceeds matrix size, choose a lower value!")

     ## And now check if everything is right...
    if( !checkData(Matrix, verbose=interactive()) )
        stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

    if( (sum(is.na(Matrix)) == 0) && (allVariables == FALSE) )
        stop("No missing values. Maybe you want to set allVariables = TRUE. Exiting\n")

    iteration = 0
    for(nPcs in evalPcs) {
       iteration = iteration + 1
       if (method == "nlpca") {
           estimate <- fitted(pca(Matrix, nPcs = nPcs, verbose = FALSE,
                                  method = method, center = TRUE,...), Matrix, nPcs = nPcs)
       } else {
           estimate <- fitted(pca(Matrix, nPcs = nPcs, verbose = FALSE,
                                  method = method, center = TRUE,...), nPcs = nPcs)
       }

        if (em == "q2") {
            # The Q2 distance
            q2 <- 1 - sum((Matrix[!missing] - estimate[!missing])^2) / sum(Matrix[!missing]^2)
            error[iteration] <- q2
        } else {
            nrmsep <- 0
            for(i in 1:ncol(Matrix)) {
                nrmsep <- nrmsep + (
                    sum((Matrix[!missing[,i], i] - estimate[!missing[,i], i])^2) /
                    (sum(!missing[,i]) * var(Matrix[,i], na.rm = TRUE))
                )
            }
            nrmsep <- nrmsep / sum(apply(missing, 2, sum) > 0)
            error[iteration] <- nrmsep
        }
        if(verbose)
           cat("The", em, "for", evalPcs[iteration], "components is:", error[iteration], "\n")
    }
            
    ret <- list()
    ret$mink <- evalPcs[which(error == min(error))]
    ret$eError <- error
    ret$evalPcs <- evalPcs

    return(ret)
}

