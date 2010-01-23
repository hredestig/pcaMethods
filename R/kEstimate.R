##' Perform cross validation to estimate the optimal number of
##' components for missing value estimation. Cross validation is
##' done for the complete subset of a variable.
##'
##' The assumption hereby is that variables that are highly correlated
##' in a distinct region (here the non-missing observations) are also
##' correlated in another (here the missing observations).  This also
##' implies that the complete subset must be large enough to be
##' representative.  For each incomplete variable, the available
##' values are divided into a user defined number of cv-segments. The
##' segments have equal size, but are chosen from a random equal
##' distribution. The non-missing values of the variable are covered
##' completely.  PPCA, BPCA, SVDimpute, Nipals PCA, llsImpute an NLPCA
##' may be used for imputation.
##'
##' The whole cross validation is repeated several times so, depending
##' on the parameters, the calculations can take very long time.  As
##' error measure the NRMSEP (see Feten et. al, 2005) or the Q2
##' distance is used.  The NRMSEP basically normalises the RMSD
##' between original data and estimate by the variable-wise
##' variance. The reason for this is that a higher variance will
##' generally lead to a higher estimation error.  If the number of
##' samples is small, the variable - wise variance may become an
##' unstable criterion and the Q2 distance should be used
##' instead. Also if variance normalisation was applied previously.
##' 
##' The method proceeds variable - wise, the NRMSEP / Q2 distance is
##' calculated for each incomplete variable and averaged
##' afterwards. This allows to easily see for wich set of variables
##' missing value imputation makes senes and for wich set no
##' imputation or something like mean-imputation should be used.  Use
##' \code{kEstimateFast} or \code{Q2} if you are not interested in
##' variable wise CV performance estimates.
##' 
##' Run time may be very high on large data sets. Especially when used
##' with complex methods like BPCA or Nipals PCA.  For PPCA, BPCA,
##' Nipals PCA and NLPCA the estimation method is called
##' \eqn{(v_{miss} \cdot segs \cdot nruncv \cdot)}{(v\_miss * segs *
##' nruncv)} times as the error for all numbers of principal
##' components can be calculated at once.  For LLSimpute and SVDimpute
##' this is not possible, and the method is called \eqn{(v_{miss}
##' \cdot segs \cdot nruncv \cdot length(evalPcs))}{(v\_miss * segs *
##' nruncv * length(evalPcs))} times. This should still be fast for
##' LLSimpute because the method allows to choose to only do the
##' estimation for one particular variable.  This saves a lot of
##' iterations.  Here, \eqn{v_{miss}}{v\_miss} is the number of
##' variables showing missing values.
##'
##' As cross validation is done variable-wise, in this function Q2 is
##' defined on single variables, not on the entire data set. This is
##' Q2 is calculated as as \eqn{\frac{\sum(x -
##' xe)^2}{\sum(x^2)}}{sum(x - xe)^2 \ sum(x^2)}, where x is the
##' currently used variable and xe it's estimate. The values are then
##' averaged over all variables. The NRMSEP is already defined
##' variable-wise. For a single variable it is then
##' \eqn{\sqrt(\frac{\sum(x - xe)^2}{(n \cdot var(x))})}{sqrt(sum(x -
##' xe)^2 \ (n * var(x)))}, where x is the variable and xe it's
##' estimate, n is the length of x.  The variable wise estimation
##' errors are returned in parameter variableWiseError.
##' @title Estimate best number of Components for missing value
##' estimation
##' @param Matrix \code{matrix} -- numeric matrix containing
##' observations in rows and variables in columns
##' @param method \code{character} -- of the methods found with
##' pcaMethods() The option llsImputeAll calls llsImpute with the
##' allVariables = TRUE parameter.
##' @param evalPcs \code{numeric} -- The principal components to use
##' for cross validation or the number of neighbour variables if used
##' with llsImpute.  Should be an array containing integer values,
##' eg. \code{evalPcs = 1:10} or \code{evalPcs = c(2,5,8)}. The NRMSEP
##' or Q2 is calculated for each component.
##' @param segs \code{numeric} -- number of segments for cross validation
##' @param nruncv \code{numeric} -- Times the whole cross validation
##' is repeated
##' @param em \code{character} -- The error measure. This can be nrmsep or q2
##' @param allVariables \code{boolean} -- If TRUE, the NRMSEP is
##' calculated for all variables, If FALSE, only the incomplete ones
##' are included. You maybe want to do this to compare several methods
##' on a  complete data set.
##' @param verbose \code{boolean} -- If TRUE, some output like the
##' variable indexes are printed to the console each iteration.
##' @param ...  Further arguments to \code{pca} or \code{nni}
##' @return   A list with:
##' \item{bestNPcs}{number of PCs or k for which the minimal average
##' NRMSEP or the maximal Q2 was obtained.}
##' \item{eError}{an array of of size length(evalPcs). Contains the
##' average error of the cross validation runs for each number of
##' components.}
##' \item{variableWiseError}{Matrix of size
##' \code{incomplete_variables} x length(evalPcs).  Contains the
##' NRMSEP or Q2 distance for each variable and each number of PCs.
##' This allows to easily see for wich variables imputation makes
##' sense and for which one it should not be done or mean imputation
##' should be used.}
##' \item{evalPcs}{The evaluated numbers of components or number of
##' neighbours  (the same as the evalPcs input parameter).}
##' \item{variableIx}{Index of the incomplete variables. This can be
##' used to map  the variable wise error to the original data.}
##' @seealso \code{\link{kEstimateFast}, \link{Q2}, \link{pca}, \link{nni}}.
##' @examples
##' ## Load a sample metabolite dataset with 5\% missing values (metaboliteData)
##' data(metaboliteData)
##' # Do cross validation with ppca for component 2:4
##' esti <- kEstimate(metaboliteData, method = "ppca", evalPcs = 2:4, nruncv=1, em="nrmsep")
##' # Plot the average NRMSEP
##' barplot(drop(esti$eError), xlab = "Components",ylab = "NRMSEP (1 iterations)")
##' # The best result was obtained for this number of PCs:
##' print(esti$bestNPcs)
##' # Now have a look at the variable wise estimation error
##' barplot(drop(esti$variableWiseError[, which(esti$evalPcs == esti$bestNPcs)]), 
##'         xlab = "Incomplete variable Index", ylab = "NRMSEP")
##' @keywords multivariate
##' @export 
##' @author Wolfram Stacklies
kEstimate <- function(Matrix, method="ppca", evalPcs=1:3, segs=3, nruncv=5,
                      em="q2", allVariables=FALSE,
                      verbose=interactive(), ...) {

  fastKE <- FALSE
  if (method == "ppca" | method == "bpca" | method == "nipals" |
      method == "nlpca")
    fastKE <- TRUE

  method <- match.arg(method, listPcaMethods())
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
    ## If the estimated observations are just scores %*% t(loadings)
    ## we can calculate all we need at once, this saves many
    ## iterations...
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

        ## We assume uniformly distributed missing values when
        ## choosing the segments
        tt <- gl(segments, ceiling(nObs / segments))[1:nObs]
        cvsegs <- split(sample(nObs), tt)
        set <- Matrix[compObs,]
        if (fastKE) {
          nrmsep <- array(0, length(evalPcs))
          q2 <- array(0, length(evalPcs))
        } else {
          nrmsep <- 0; q2 <- 0
        }
        
        for (i in 1:length(cvsegs)) {
          n <- length(cvsegs[[i]]) # n is the number of created
                                   # missing values
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
              ## Error of prediction, error is calculated for removed
              ## elements only
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
            ## Error of prediction, error is calculated for removed
            ## elements only
            if (em == "nrmsep") {
              nrmsep <- nrmsep + sum( (original - estimate)^2)
            } else {
              q2 <- q2 + sum( (original - estimate)^2 )
            }
          }
        } ## iteration over cv segments
        
        if (fastKE) {
          if (em == "nrmsep") {
            error[pos, ] <-
              error[pos, ] + nrmsep / (nrow(set) * var(set[,index]))
          } else
          error[pos, ] <- error[pos, ] + (1 - (q2 / sum(set[, index]^2)))
        } else {
          if (em == "nrmsep") {
            error[pos, iteration] <- error[pos, iteration] + 
              nrmsep / (nrow(set) * var(set[,index]))
          } else
          error[pos, iteration] <-
            error[pos, iteration] + (1 - (q2 / sum(set[, index]^2)))
        }
      } # iteration over variables
      if (verbose) cat("\n")
      
    } #iteration over nruncv
    
    ## The error is the sum over the independent cross validation runs
    error <- error / nruncv
    
    if (verbose && !fastKE)
      cat("The average", em, "for k =", iteration, "is", 
          sum(error[,iteration]) / nrow(error), "\n")

    ## if nlpca, ppca, bpca, nipals we do not need to iterate over the
    ## number of components...
    if (fastKE) break
  } # iteration over number components

  if (em == "nrmsep")
    avgError <- sqrt(apply(error, 2, sum) / nrow(error))
  else
    avgError <- apply(error, 2, sum) / nrow(error)

  ret <- list()
  if (em == "nrmsep")
    ret$bestNPcs <- evalPcs[which(avgError == min(avgError))]
  else ret$bestNPcs <- evalPcs[which(avgError == max(avgError))]
  ret$eError <- avgError
  if(em == "nrmsep") ret$variableWiseError <- sqrt(error)
  else ret$variableWiseError <- error
  ret$evalPcs <- evalPcs
  ret$variableIx <- missIx

  return(ret)
}
