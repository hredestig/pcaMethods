##' ONB = orth(mat) is an orthonormal basis for the range of matrix
##' mat.  That is, ONB' * ONB = I, the columns of ONB span the same
##' space as the columns of mat, and the number of columns of ONB is
##' the rank of mat.
##' @title Calculate an orthonormal basis
##' @param mat matrix to calculate orthonormal base
##' @param skipInac  do not include components with precision below
##' .Machine$double.eps if TRUE
##' @return orthonormal basis for the range of matrix
##' @author Wolfram Stacklies
orth <- function(mat, skipInac = FALSE) {

  if(nrow(mat) > ncol(mat)) {
    leftSVs <- ncol(mat)
  } else {
    leftSVs <- nrow(mat)
  }

  result <- svd(mat, nu =  leftSVs, nv = ncol(mat))
  U <- result[[2]]
  S <- result[[1]]
  V <- result[[3]]

  m <- nrow(mat)
  n <- ncol(mat)

  if(m > 1) { 
    s <- diag(S, nrow = length(S))
  } else     if(m == 1) { 
    s <- S[1] 
  } else { 
    s <- 0 
  }

  tol <- max(m,n) * max(s) * .Machine$double.eps
  r <- sum(s > tol)
  if ( r < ncol(U) ) {
    if (skipInac) {
      warning("Precision for components ", r + 1 , " - ", ncol(U), 
              " is below .Machine$double.eps. \n",
              "Results for those components are likely to be inaccurate!!\n",
              "These component(s) are not included in the returned solution!!\n")
    } else {
      warning("Precision for components ", r + 1 , " - ", ncol(U), 
              " is below .Machine$double.eps. \n",
              "Results for those components are likely to be inaccurate!!\n")
    }
  }
  
  if (skipInac) {
    ONB <- U[, 1:r, drop=FALSE]
    ## Assing correct row and colnames
    rownames(ONB) <- labels(mat[, 1:r, drop=FALSE])[[1]];
    colnames(ONB) <- labels(mat[, 1:r, drop=FALSE])[[2]];
  } else {
    ONB<-U
    ## Assing correct row and colnames
    rownames(ONB) <- labels(mat)[[1]];
    colnames(ONB) <- labels(mat)[[2]];
  }

  return(ONB)
}
