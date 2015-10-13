##' Creates a large matrix B consisting of an M-by-N tiling of copies
##' of A
##' @title Replicate and tile an array.
##' @param mat numeric matrix
##' @param M number of copies in vertical direction
##' @param N number of copies in horizontal direction
##' @return Matrix consiting of M-by-N tiling copies of input matrix
##' @author Wolfram Stacklies
repmat <- function(mat, M, N) {

  ## Check if all input parameters are correct
  if( !all(M > 0, N > 0) ) {
    stop("M and N must be > 0")
  }    
  
  ## Convert array to matrix
  ma <- mat
  if(!is.matrix(mat)) {
    ma <- matrix(mat, nrow=1)
  }

  rows <- nrow(ma)
  cols <- ncol(ma)
  replicate <- matrix(0, rows * M, cols * N)

  for (i in 1:M) {
    for(j in 1:N) {
      start_row <- (i - 1) * rows + 1
      end_row <- i * rows
      start_col <- (j - 1) * cols + 1
      end_col <- j * cols
      replicate[start_row:end_row, start_col:end_col] <- ma
    }
  }

  return(replicate)
}

