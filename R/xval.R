##' Internal cross-validation can be used for estimating the level of
##' structure in a data set and to optimise the choice of number of
##' principal components.
##'
##' This method calculates \eqn{Q^2} for a PCA model. This is the
##' cross-validated version of \eqn{R^2} and can be interpreted as the
##' ratio of variance that can be predicted independently by the PCA
##' model. Poor (low) \eqn{Q^2} indicates that the PCA model only
##' describes noise and that the model is unrelated to the true data
##' structure. The definition of \eqn{Q^2} is: \deqn{Q^2=1 -
##' \frac{\sum_{i}^{k}\sum_{j}^{n}(x -
##' \hat{x})^2}{\sum_{i}^{k}\sum_{j}^{n}x^2}}{Q^2=1 - sum_i^k
##' sum_j^n (x - \hat{x})^2 / \sum_i^k \sum_j^n(x^2)} for the matrix
##' \eqn{x} which has \eqn{n} rows and \eqn{k} columns. For a given
##' number of PC's x is estimated as \eqn{\hat{x}=TP'} (T are scores
##' and P are loadings). Although this defines the leave-one-out
##' cross-validation this is  not what is performed if fold is less
##' than the number of rows and/or columns.  In 'impute' type CV,
##' diagonal rows of elements in the matrix are deleted and the
##' re-estimated.  In 'krzanowski' type CV, rows are sequentially left
##' out to build fold PCA models which give the loadings. Then,
##' columns are sequentially left out to build fold models for
##' scores. By combining scores and loadings from different models, we
##' can estimate completely left out values.  The two types may seem
##' similar but can give very different results, krzanowski typically
##' yields more stable and reliable result for estimating data
##' structure whereas impute is better for evaluating missing value
##' imputation performance. Note that since Krzanowski CV operates on
##' a reduced matrix, it is not possible estimate Q2 for all
##' components and the result vector may therefore be shorter than
##' \code{nPcs(object)}.
##' @title Cross-validation for PCA
##' @param object A \code{pcaRes} object (result from previous PCA
##' analysis.)
##' @param originalData The matrix (or ExpressionSet) that used to
##' obtain the pcaRes object. 
##' @param fold The number of groups to divide the data in.
##' @param nruncv The number of times to repeat the whole
##' cross-validation
##' @param type krzanowski or imputation type cross-validation
##' @param verbose \code{boolean} If TRUE Q2 outputs a primitive
##' progress bar.
##' @param variables indices of the variables to use during
##' cross-validation calculation. Other variables are kept as they are
##' and do not contribute to the total sum-of-squares.
##' @param ... Further arguments passed to the \code{\link{pca}} function called
##' within Q2.
##' @return A matrix or vector with \eqn{Q^2} estimates.
##' @export
##' @references Krzanowski, WJ. Cross-validation in principal
##' component analysis. Biometrics. 1987(43):3,575-584
##' @examples
##' data(iris)
##' x <- iris[,1:4]
##' pcIr <- pca(x, nPcs=3)
##' q2 <- Q2(pcIr, x)
##' barplot(q2, main="Krzanowski CV", xlab="Number of PCs", ylab=expression(Q^2))
##' ## q2 for a single variable
##' Q2(pcIr, x, variables=2)
##' pcIr <- pca(x, nPcs=3, method="nipals")
##' q2 <- Q2(pcIr, x, type="impute")
##' barplot(q2, main="Imputation CV", xlab="Number of PCs", ylab=expression(Q^2))
##' @author Henning Redestig, Ondrej Mikula
##' @keywords multivariate
Q2 <- function (object, originalData=completeObs(object), fold=5, 
                nruncv=1, type=c("krzanowski", "impute"), verbose=interactive(),
                variables=1:nVar(object), ...) {
  type <- match.arg(type)
  if (inherits(originalData, "ExpressionSet")) {
    set <- originalData
    originalData <- t(exprs(originalData))
  }
  if (is.null(originalData)) 
    stop("missing data when estimating Q2")
  originalData <- as.matrix(originalData)
  originalData <- prep(originalData, scale=scl(object), center=center(object))
  nR <- nObs(object)
  nC <- nVar(object)
  if (nR != nrow(originalData) | nC != ncol(originalData)) 
    stop("data and model dimensions do not match")
  if (fold > max(nR, nC)) 
    stop("fold must be equal or less to max dimension of original data")
  if (method(object) %in% c("svd") & type != "krzanowski") 
    stop("Chosen PCA method must use krzanowski type cv")
  if (method(object) %in% c("llsImpute") & type != "impute") 
    stop("Chosen PCA method must use impute type cv")
  
  if (is.logical(variables))
    variables <- which(variables)

  ssx <- sum(originalData[, variables]^2, na.rm=TRUE)
  for (nr in 1:nruncv) {
    if (type == "impute") {
      nP <- nPcs(object)
      press <- rep(0, nP)
      q2 <- matrix(NA, nP, ncol=nruncv)
      seg <- list()
      nDiag <- max(nR, nC)
      diagPerFold <- floor(nDiag/fold)
      suppressWarnings(diags <- matrix(1:nDiag, nrow=diagPerFold, ncol=fold, byrow=TRUE))
      if (diagPerFold == 0 || diagPerFold > (nDiag/2)) 
        stop("Matrix could not be safely divided into ", 
             fold, " segments. Choose a different fold or provide the desired segments")
      if (nDiag%%fold > 0) 
        warning("Validation incomplete: ", (nDiag%%fold) * 
                  min(dim(originalData)),
                " values were left out of from cross validation, Q2 estimate will be biased.")
      for (i in 1:ncol(diags))
        seg[[i]] <- which(is.na(deletediagonals(originalData, diags[, i])))
      if (verbose) {
        message("Doing ", length(seg), " fold ", "cross validation")
        pb <- txtProgressBar(0, length(seg), style=3, width=20)
      }
      j <- 0
      for (i in seg) {
        j <- j + 1
        if (verbose) 
          setTxtProgressBar(pb, j)
        test <- originalData
        test[i] <- NA
        test <- tempFixNas(test)
        if (method(object) != "llsImpute") {
          pc <- pca(test, nPcs=nP, method=method(object), 
					verbose=FALSE, center=centered(object), scale=object@scaled, ...)
        }
        for (np in 1:nP) {
          if (method(object) == "llsImpute") {
            fittedData <- completeObs(llsImpute(test, k=np, allVariables=TRUE, center=FALSE))
          }
          else {
            if (method(object) == "nlpca") 
              fittedData <- fitted(pc, data=test, nPcs=np)
            else fittedData <- fitted(pc, data=NULL, nPcs=np)
          }
          ii <- i[ceiling(i / nR) %in% variables]
          press[np] <- press[np] + sum((originalData[ii] - fittedData[ii])^2, na.rm=TRUE)
        }
      }
    }

    if (type == "krzanowski") {
      rseg <- split(sample(1:nR), rep(1:fold, ceiling(nR/fold))[1:nR])
      cseg <- split(sample(1:nC), rep(1:fold, ceiling(nC/fold))[1:nC])
      nP <- min(nR - max(sapply(rseg, length)), nC - max(sapply(cseg, length)), nPcs(object))
      q2 <- matrix(NA, nP, ncol=nruncv)
      press <- rep(0, nP)
      foldC <- length(cseg)
      foldR <- length(rseg)
      tcv <- array(0, dim=c(foldC, nR, nP))
      pcv <- array(0, dim=c(foldR, nC, nP))
      for (f in 1:foldC) {
        test <- tempFixNas(originalData[, -cseg[[f]]])
        tcv[f, , ] <- scores(pca(test, nPcs=nP, method=method(object), 
                                 verbose=FALSE, center=centered(object), scale=object@scaled, ...))
        for (p in 1:nP) {
          if (cor(tcv[f, , p], scores(object)[, p]) < 0) 
            tcv[f, , p] <- tcv[f, , p] * -1
        }
      }
      for (f in 1:foldR) {
        test <- tempFixNas(originalData[-rseg[[f]], ])
        pcv[f, , ] <- loadings(pca(test, nPcs=nP, method=method(object), 
                                   verbose=FALSE, center=centered(object), scale=object@scaled, ...))
        for (p in 1:nP) {
          if (cor(pcv[f, , p], loadings(object)[, p]) < 0) 
            pcv[f, , p] <- pcv[f, , p] * -1
        }
      }
      press <- rep(0, nP)
      for (p in 1:nP)
        for (fr in 1:foldR)
          for (fc in 1:foldC) 
            press[p] <- press[p] + sum((
              originalData[rseg[[fr]], cseg[[fc]]] - 
                (tcv[fc, , ][, 1:p, drop=FALSE] %*% t(pcv[fr, , ][,1:p, drop=FALSE]))
              [rseg[[fr]], intersect(cseg[[fc]], variables)])^2, na.rm=TRUE)
    }
    q2[, nr] <- 1 - press/ssx
  }

  if (verbose) 
    message("\n")
  rownames(q2) <- paste("PC", 1:nrow(q2))
  drop(q2)
}

##' Simply replace completely missing rows or cols with zeroes.
##' @title Temporary fix for missing values
##' @param mat a matrix
##' @return The original matrix with completely missing rows/cols
##' filled with zeroes.
##' @author Henning Redestig 
tempFixNas <- function(mat) {
  badRows <- apply(mat, 1, function(x) all(is.na(x)))
  badCols <- apply(mat, 2, function(x) all(is.na(x)))
  mat[ badRows,] <- 0
  mat[,badCols ] <- 0
  mat
}

##' Replace a diagonal of elements of a matrix with NA
##'
##' Used for creating artifical missing values in matrices without
##' causing any full row or column to be completely missing
##' @title Delete diagonals
##' @param x The matrix  
##' @param diagonals The diagonal to be replaced, i.e. the first,
##' second and so on when looking at the fat version of the matrix
##' (transposed or not) counting from the bottom.
##' Can be a vector to delete more than one diagonal.
##' @return The original matrix with some values missing
##' @author Henning Redestig 
deletediagonals <- function(x, diagonals=1) {
  wastransposed <- FALSE
  if (dim(x)[1] > dim(x)[2]) {          # matrix must be lying down
    x <- t(x)
    wastransposed <- TRUE
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (!all(diagonals <= nc)) {
    stop(paste("Order of diagonal number", max(diagonals),  "is out of bound"))
  }
  indexmatrix <- matrix(1 : (nr * nc), ncol=nc, nrow=nr)
  finalmatrix <- matrix(ncol=(nr - 1 + nc), nrow=nr)
  finalmatrix[,1 : (nr - 1)] <- indexmatrix[,rev((nc : 1)[1 : (nr - 1)])]
  finalmatrix[,nr : (nr - 1 + nc)] <- indexmatrix
  dia <- 1 + 0:(nr - 1) * (nr + 1)
  finalIndices <- NULL
  for (i in 1:length(diagonals)) {
    indicestodelete <- finalmatrix[dia + (diagonals[i] - 1) * nr]
    x[indicestodelete] <- NA
    finalIndices <- c(finalIndices, indicestodelete)
  }
  if (wastransposed) x <- t(x)
  return(x)
}

##' Get cross-validation segments that have (as far as possible) the
##' same ratio of all classes (if classes are present)
##' @title Get CV segments
##' @param x a factor, character or numeric vector that describes
##' class membership of a set of items, or, a numeric vector
##' indicating unique indices of items, or, a numeric of length 1 that
##' describes the number of items to segment (without any classes)
##' @param fold the desired number of segments
##' @param seed randomization seed for reproducibility
##' @return a list where each element is a set of indices that defines
##' the CV segment.
##' @examples
##' seg <- cvseg(iris$Species, 10)
##' sapply(seg, function(s) table(iris$Species[s]))
##' cvseg(20, 10)
##' @seealso the \code{cvsegments} function in the \code{pls} package
##' @export
##' @author Henning Redestig
cvseg <- function(x, fold=7, seed=NULL) {
  if(any(table(x) > 1)) {
    if(any(table(x) < fold)) {
      fold <- min(table(x))
    }
    if(fold < 2)
      stop("too few observations in the smallest class")
    res <-
      sapply(unique(x), function(z) {
        if(!is.null(seed))
          set.seed(seed)
        tmp <- sample(which(x == z))
        seg <- matrix(c(tmp,
                        rep(NA, ifelse(length(tmp) %% fold ==0, 0,
                                       fold - (length(tmp) %% fold)))),
                      nrow=fold)
      },simplify=FALSE)
    res <- do.call("cbind", res)
  } else {
    if(length(x) == 1) x <- 1:x
    res <- matrix(sample(c(x,
                           rep(NA, ifelse(length(x) %% fold ==0, 0,
                                          fold - (length(x) %% fold))))),
                  nrow=fold)
  }
  res <- res[!apply(is.na(res), 1, all),,drop=FALSE]
  res
  lapply(as.data.frame(t(res)), function(x) c(na.omit(x)))
}
