#####################################################################################
## repmat <- function(mat, M, N)
##
## Replicate and tile an array.
##    B = repmat(A,M,N) creates a large matrix B consisting of an M-by-N
##    tiling of copies of A
##
## Parameters:
## matrix   - numeric matrix
## M        - number of copies in vertical direction
## N        - number of copies in horizontal direction
##
## Return:
## replicate    - Matrix consiting of M-by-N tiling copies of input matrix
##
## Author: Wolfram Stacklies
##         Max Planck Institut fuer Molekulare Pflanzenphysiologie
##         Golm, Germany
## Date:   11.04.2006
##
## Contact:    wolfram.stacklies@gmail.com
##
#####################################################################################

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

