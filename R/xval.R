Q2 <- function(object, originalData, nPcs=object@nPcs, fold=5, nruncv=10) {

  originalData <- as.matrix(originalData)
  
  if(length(object@subset) > 0)
    originalData <- originalData[,object@subset]
  
  originalData <- prep(originalData, center=object@centered, scale=object@scaled)

  nR <- nrow(originalData)
  nC <- ncol(originalData)

  ssx <- sum(originalData^2, na.rm=TRUE)
  
  if(object@method %in% c("svd"))
    stop("Chosen PCA method can not handle missing values!")
  result <- matrix(NA, nPcs, ncol=nruncv)
  for(nr in 1:nruncv) {
    ## Error here, no complete rows and no complete columns may be deleted
    seg <- cvsegments(length(originalData), fold)

    ## Choose a segment such that no rows or columns are entirely filled with NA values
    for (i in 1:fold) {
      test <- originalData
      test[ seg[[i]] ] <- NA
      ## Check if we created rows or cols wich are entirely NA
      count = 0
      while ( (sum( apply(is.na(test), 1, sum) == ncol(test) ) > 0) ||
             (sum( apply(is.na(test), 2, sum) == nrow(test) ) > 0) ) {
        ## Now choose a new segment
        segNew <- cvsegments(length(originalData), fold)
        seg[i] = segNew[i]
        test <- originalData
        test[ seg[[i]] ] <- NA
        count = count + 1
        if(count > 100) {
          stop("cvsegments: No valid set found during 100 iterations, exiting. 
					Choose a lower fold value!!")
        }
      }
    }
    
    press <- rep(0, nPcs)
    for(i in seg) {
      test <- originalData
      test[i] <- NA
      pc <- pca(test, nPcs=nPcs, method=object@method, verbose=FALSE,
                center=object@centered, scale=object@scaled)
      T <- pc@scores
      P <- pc@loadings
      
      ## add up to the press estimate
      for(np in 1:nPcs) 
        press[np] <- press[np] +
          sum((originalData[i] -
               (T[,1:np, drop=FALSE] %*% t(P[,1:np, drop=FALSE]))[i])^2, na.rm=TRUE)
    }
    ## now calculate the q2 estimate for each pc
    result[,nr] <- 1 - press / ssx
  }
  rownames(result) <- paste("PC", 1:nrow(result))
  result
}
