Q2 <- function(object, originalData, nPcs=object@nPcs, fold=5, nruncv=10, segments=NULL, verbose=interactive(),...) {

  originalData <- as.matrix(originalData)
  
  if(length(object@subset) > 0)
    originalData <- originalData[,object@subset]
  
  originalData <- prep(originalData, center=object@centered, scale=object@scaled)

  nR <- nrow(originalData)
  nC <- ncol(originalData)

  ssx <- sum(originalData^2, na.rm=TRUE)
  
  if(object@method %in% c("svd"))
    stop("Chosen PCA method can not handle missing values and missing value free cross validation is not supported")
  result <- matrix(NA, nPcs, ncol=nruncv)
  for(nr in 1:nruncv) {
    seg <- segments
    if(is.null(seg)) {

      nDiag <- max(dim(originalData))
      diagPerFold <- floor(nDiag / fold)
      suppressWarnings(diags <- matrix(1:nDiag, nrow=diagPerFold, ncol=fold, byrow=TRUE))
      if(diagPerFold == 0 || diagPerFold > (nDiag / 2))
        stop("Matrix could not be safely divided into ", fold, " segment(s). Choose a different fold or provide the desired segments")
      if(nDiag %% fold > 0)
        warning("Validation incomplete: ", (nDiag %% fold) * min(dim(originalData)), " value(s) were left out of from cross validation, Q2 estimate will be biased.")

      ## not get the indices of those diagonals
      for(i in 1:ncol(diags)) 
        seg[[i]] <- which(is.na(deletediagonals(originalData, diags[,i])))

      ## <to remove later>
      ## check that we didnt make a mistake
      for(i in 1:length(seg)) {
        tmp <- originalData
        tmp[seg[[i]]] <- NA
        if(any(apply(tmp, 1, function(x) sum(is.na(x))) == ncol(tmp)) ||
           any(apply(tmp, 2, function(x) sum(is.na(x))) == nrow(tmp)))
          stop("ooops! a column or a row was completely lost, this should not have happened. Please contact the maintainer")
      }
      ## </to remove later>
      
    }
    
    press <- rep(0, nPcs)
      if(verbose)
        message("Doing ", length(seg), " fold ", "cross validation")
    for(i in seg) {
      if(verbose)
        cat(".")
      test <- originalData
      test[i] <- NA
      pc <- pca(test, nPcs=nPcs, method=object@method, verbose=FALSE,
                center=object@centered, scale=object@scaled,...)
      
      ## add up to the press estimate
      for(np in 1:nPcs) {
        if(object@method == "nlpca")
          fittedData <- fitted(pc, data=test, nPcs=np)
        else
          fittedData <- fitted(pc, data=NULL, nPcs=np)
        press[np] <- press[np] +
          sum((originalData[i] - fittedData[i])^2, na.rm=TRUE)
      }
    }
    if(verbose)
      cat("\nDone\n")
    ## now calculate the q2 estimate for each pc
    result[,nr] <- 1 - press / ssx
  }
  rownames(result) <- paste("PC", 1:nrow(result))
  result
}

deletediagonals <- function(x, diagonals=1) {

  ##<..Beg Rdocu..>
  ## ~name~
  ##   deletediagonals
  ## ~title~
  ##   Replace a diagonal of elements of a matrix with NA
  ## ~description~
  ##   Used for creating artifical missing values in matrices without
  ##   causing any full row or column to be completely missing
  ## ~usage~
  ##   deletediagonals(x, diagonals=1)
  ## ~arguments~
  ##   ~-x~
  ##     The matrix
  ##   ~-diagonals~
  ##     The diagonal to be replaced, i.e. the first, second and so on
  ##     when looking at the fat version of the matrix (transposed or
  ##     not) counting from the bottom. Can be a vector to delete more
  ##     than one diagonal.
  ## ~details~
  ##   
  ## ~value~
  ##   The original matrix with some values missing
  ## ~seealso~
  ##   cvsegments from pls
  ## ~examples~
  ##   deletediagonals(iris[,1:4], 1)
  ## ~keywords~
  ##   
  ## ~author~
  ##   Henning Redestig <redestig[at]mpimp-golm.mpg.de>
  ##>..End Rdocu..<
  
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
  indexmatrix <- matrix(1 : (nr * nc), ncol = nc, nrow = nr)
  finalmatrix <- matrix(ncol = (nr - 1 + nc), nrow = nr)
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
