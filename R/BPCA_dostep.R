##' The function contains the actual implementation of the BPCA
##' component  estimation. It performs one step of the BPCA EM
##' algorithm.  It is called 'maxStep' times from within the main loop
##' in BPCAestimate.
##'
##' This function is NOT intended to be run standalone.
##' @title Do BPCA estimation step
##' @param M Data structure containing all needed information. See the
##' source documentation of BPCA_initmodel for details
##' @param y Numeric original data matrix
##' @return Updated version of the data structure
##' @author Wolfram Stacklies
BPCA_dostep <- function(M,y) {

  ## Empty matrix in which the scores are copied
  M$scores <- matrix(NA, M$rows, M$comps)

  ## Expectation step for data without missing values
  Rx <- diag(M$comps) + M$tau * t(M$PA) %*% M$PA + M$SigW
  Rxinv <- solve(Rx)
  idx <- M$row_nomiss

  if (length(idx) == 0) {
    trS <- 0
    T <- 0
  } else {
    dy <- y[idx,, drop=FALSE] - repmat(M$mean, length(idx), 1)
    x <- M$tau * Rxinv %*% t(M$PA) %*% t(dy)
    T <- t(dy) %*% t(x)
    trS <- sum(sum(dy * dy))

    ## Assign the scores for complete rows
    xTranspose <- t(x)
    for (i in 1:length(idx)) {
      M$scores[idx[i],] <- xTranspose[i,]
    }
  }
  ## Expectation step for incomplete data
  if( length(M$row_miss) > 0) {
    for(n in 1:length(M$row_miss)) {
      i  <- M$row_miss[n]
      dyo <- y[ i, !M$nans[i,], drop=FALSE] - M$mean[ !M$nans[i,], drop=FALSE]
      Wm <- M$PA[ M$nans[i,],, drop=FALSE]
      Wo <- M$PA[ !M$nans[i,],, drop=FALSE]
      Rxinv <- solve( (Rx - M$tau * t(Wm) %*% Wm))
      ex  <- M$tau * t(Wo) %*% t(dyo)
      x <- Rxinv %*% ex
      dym <- Wm %*% x
      dy <- y[i,, drop=FALSE]
      dy[ !M$nans[i,] ] <- t(dyo)
      dy[ M$nans[i,] ] <- t(dym)
      M$yest[i,] <- dy + M$mean
      T <- T + t(dy) %*% t(x)
      T[ M$nans[i,], ] <- T[ M$nans[i,],, drop=FALSE] + Wm %*% Rxinv
      trS <- trS + dy %*% t(dy) + sum(M$nans[i,]) / M$tau + 
        sum( diag(Wm %*% Rxinv %*% t(Wm)) ) 
      trS <- trS[1,1]
      ## Assign the scores for rows containing missing values
      M$scores[M$row_miss[n],] <- t(x)
    }
  }
  T <- T / M$rows
  trS <- trS / M$rows

  ## Maximation step
  Rxinv <- solve(Rx)
  Dw <- Rxinv + M$tau * t(T) %*% M$PA %*% Rxinv + 
    diag(M$alpha, nrow = length(M$alpha)) / M$rows
  Dwinv <- solve(Dw)
  M$PA <- T %*% Dwinv ## The new estimate of the principal axes (loadings)
  M$tau <- (M$cols + 2 * M$gtau0 / M$rows) / (trS - sum(diag(t(T) %*% M$PA)) +
                                              (M$mean %*% t(M$mean) * M$gmu0 + 2 * M$gtau0 / M$btau0) / M$rows)
  M$tau <- M$tau[1,1] ## convert to scalar
  M$SigW <- Dwinv * (M$cols / M$rows)
  M$alpha <- (2 * M$galpha0 + M$cols) / (M$tau * diag(t(M$PA) %*% M$PA) + 
                                         diag(M$SigW) + 2 * M$galpha0 / M$balpha0)

  return(M)
}


