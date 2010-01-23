##' Line search for conjugate gradient
##' @param nlnet The nlnet
##' @param dw ..
##' @param e0 ..
##' @param ttGuess ..
##' @param trainIn Training data
##' @param trainOut Fitted data
##' @param verbose  logical, print messages
##' @return ...
##' @author Henning Redestig, Matthias Scholz
lineSearch <- function(nlnet, dw, e0, ttGuess, trainIn, trainOut, verbose) {
  iterGoldenSectionSearch <- 6
  alpha <- 0.618034
  tt <- rep(0, 4)
  e <- rep(0, 4)
  tmpnlnet <- forkNlpcaNet(nlnet)
  

  tt[1] <- 0
  e[1] <- e0
  tt[4] <- ttGuess
  tmpnlnet@weights$set(nlnet@weights$current() + tt[4] * dw)
  e[4] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error

  if(e[4] > e[1]) {                     #got final interval calculate tt[2] and tt[3]
    tt[2] <- tt[1] + (1 - alpha) * (tt[4] - tt[1])
    tmpnlnet@weights$set(nlnet@weights$current() + tt[2] * dw)
    e[2] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
    tt[3] <- tt[1] + alpha * (tt[4] - tt[1])
    tmpnlnet@weights$set(nlnet@weights$current() + tt[3] * dw)
    e[3] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
  }
  else {                                #expand, add new tt[4]
    tt[3] <- tt[4]
    e[3] <- e[4]
    tt[4] <- (1 + alpha) * tt[4]
    tmpnlnet@weights$set(nlnet@weights$current() + tt[4] * dw)
    e[4] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
    if(e[4] > e[3]) {                   #got final interval, calculate tt[2]
      tt[2] <- tt[1] + (1 - alpha) * (tt[4] - tt[1])
      tmpnlnet@weights$set(nlnet@weights$current() + tt[2] * dw)
      e[2] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
    }
    else {                              #expand: add new tt[4]
      i <- 1
      while(e[4] < e[3] && i < 50) {
        tt[2] <- tt[3]
        e[2] <- e[3]
        tt[3] <- tt[4]
        e[3] <- e[4]
        tt[4] <- (1 + alpha) * tt[4]
        tmpnlnet@weights$set(nlnet@weights$current() + tt[4] * dw)
        e[4] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
        i <- i + 1
        if(verbose && i == 50)
          cat("^")
      }
    }
  }

  ## golden section search

  for(i in 1:iterGoldenSectionSearch) {
    if(e[3] > e[2]) {
      tt[4] <- tt[3]                      #remove right value tt[4]
      e[4] <- e[3]
      tt[3] <- tt[2]
      e[3] <- e[2]
      tt[2] <- tt[1] + (1 - alpha) * (tt[4] - tt[1]) #split left interval
      tmpnlnet@weights$set(nlnet@weights$current() + tt[2] * dw)
      e[2] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
    }
    else {
      tt[1] <- tt[2]                      #remove left t value tt[1]
      e[1] <- e[2]
      tt[2] <- tt[3]
      e[2] <- e[3]
      tt[3] <- tt[1] + alpha * (tt[4] - tt[1]) #split right interval
      tmpnlnet@weights$set(nlnet@weights$current() + tt[3] * dw)
      e[3] <- nlnet@error(tmpnlnet, trainIn, trainOut)$error
    }
  }
  if(e[2] < e[3]) {
    eBest <- e[2]
    ttBest <- tt[2]
  }
  else {
    eBest <- e[3]
    ttBest <- tt[3]
  }

  wBest <- nlnet@weights$current() + ttBest * dw
  return(list(wBest=wBest, eBest=eBest, ttBest=ttBest))
}

##' Linear kernel
##' @param x datum
##' @return Input value
##' @author Henning Redestig, Matthias Scholz
linr <- function(x) x
