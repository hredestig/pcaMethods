##' This is a simple estimator for the optimal number of componets
##' when applying PCA or LLSimpute for missing value estimation.  No
##' cross validation is performed, instead the estimation quality is
##' defined as Matrix[!missing] - Estimate[!missing]. This will give a
##' relatively rough estimate, but the number of iterations equals the
##' length of the parameter evalPcs.\cr Does not work with LLSimpute!!
##' As error measure the NRMSEP (see Feten et. al, 2005) or the Q2
##' distance is used.  The NRMSEP basically normalises the RMSD
##' between original data and estimate by the variable-wise
##' variance. The reason for this is that a higher variance will
##' generally lead to a higher estimation error.  If the number of
##' samples is small, the gene - wise variance may become an unstable
##' criterion and the Q2 distance should be used instead. Also if
##' variance normalisation was applied previously.
##' @title Estimate best number of Components for missing value estimation
##' @param Matrix \code{matrix} -- numeric matrix containing
##' observations in rows and variables in columns
##' @param method \code{character} -- a valid pca method (see
##' \code{\link{pca}}).
##' @param evalPcs \code{numeric} -- The principal components to use
##' for cross validation or cluster sizes if used with
##' llsImpute. Should be an array containing integer values,
##' eg. evalPcs = 1:10 or evalPcs = C(2,5,8).The NRMSEP is calculated
##' for each component.
##' @param em \code{character} -- The error measure. This can be
##' nrmsep or q2
##' @param allVariables \code{boolean} -- If TRUE, the NRMSEP is
##' calculated for all variables, If FALSE, only the incomplete ones
##' are included. You maybe want to do this to compare several methods
##' on a  complete data set.
##' @param verbose \code{boolean} -- If TRUE, the NRMSEP and the
##' variance are printed to the console each iteration.
##' @param ... Further arguments to \code{pca}
##' @return   \item{list}{Returns a list with the elements:
##' \itemize{
##' \item minNPcs - number of PCs for which the minimal average NRMSEP
##' was obtained
##' \item eError - an array of of size length(evalPcs). Contains the
##' estimation error for each number of
##' components.
##' \item evalPcs - The evaluated numbers of components or
##' cluster sizes  (the same as the evalPcs input parameter). }}
##' @seealso \code{\link{kEstimate}}.
##' @export
##' @examples
##' data(metaboliteData)
##' # Estimate best number of PCs with ppca for component 2:4
##' esti <- kEstimateFast(t(metaboliteData), method = "ppca", evalPcs = 2:4, em="nrmsep")
##' barplot(drop(esti$eError), xlab = "Components",ylab = "NRMSEP (1 iterations)")
##' # The best k value is:
##' print(esti$minNPcs)
##' @keywords multivariate
##' @author Wolfram Stacklies
kEstimateFast <- function(Matrix, method = "ppca", evalPcs = 1:3,
                          em = "nrmsep", allVariables = FALSE,
                          verbose = interactive(), ...) {

    method <-
      match.arg(method, c("ppca", "bpca", "svdImpute", "nipals", "nlpca"))
    em <- match.arg(em, c("nrmsep", "q2"))
    maxPcs <- max(evalPcs)
    lengthPcs <- length(evalPcs)
    missing <- is.na(Matrix)
    error <- array(0, lengthPcs)

     ## If the data is a data frame, convert it into a matrix
    Matrix <- as.matrix(Matrix, rownames.force=TRUE)
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
    if (em == "nrmsep") ret$bestNPcs <- evalPcs[which(error == min(error))]
    else ret$bestNPcs <- evalPcs[which(error == max(error))]
    ret$eError <- error
    ret$evalPcs <- evalPcs

    return(ret)
}

