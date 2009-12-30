## Hope the fixes a strange bug during package test on bioconductor...
.onLoad <- function(libname, pkgname) {
  require("methods")
}

pca <- function(object, method=c("svd", "nipals", "bpca", "ppca",
                          "svdImpute", "nlpca", "robustPca"),
                subset=numeric(),...) {
  
  isExprSet <- FALSE
  if(inherits(object, "ExpressionSet")) {
    set <- object
    isExprSet <- TRUE
    object <- t(exprs(object))
  }

  method <- match.arg(method)

  ## Do some basic checks of the data. We exit if the data contains
  ## NaN or Inf values or is not numeric.
  if ( !checkData(as.matrix(object), verbose=interactive()) )
    stop("Invalid data format, exiting...\n",
         "Run checkData(data, verbose=TRUE) for details\n")
  
  missing <- sum(is.na(object))

  if(length(subset) > 0)
    object <- object[,subset]

  if(missing > 0 & method == "svd") {
    warning("Found missing values, using the nipals method instead of requested svd")
    method <- "nipals"
  }
  
  switch(method,
         svd = {
           res <- svdPca(object, ...)
         },
         nipals = {
           res <- nipalsPca(object, ...)
         },
         bpca = {
           res <- bpca(as.matrix(object),...)
         },
         ppca = {
           res <- ppca(as.matrix(object),...)
         },
         svdImpute = {
           res <- svdImpute(as.matrix(object),...)
         },
         robustPca = {
           res <- robustPca(as.matrix(object), ...)
         },
         nlpca = {
           res <- nlpca(as.matrix(object),...)
         })

  return(res)
}


##
## This method currently only serves for llsImpute
##
nni <- function(object, method=c("llsImpute"), subset=numeric(), ...) {

  isExprSet <- FALSE
  if(inherits(object, "ExpressionSet")) {
    set <- object
    isExprSet <- TRUE
    object <- t(exprs(object))
  }

  method <- match.arg(method)

                                        ## Do some basic checks of the data. We exit if the data contains NaN or Inf
                                        ## values or is not numeric.
  if ( !checkData(as.matrix(object), verbose = interactive()) )
    stop("Invalid data format, exiting...\n",
         "Run checkData(data, verbose = TRUE) for details\n")

  missing <- sum(is.na(object))
  if(length(subset) > 0)
    object <- object[,subset]

  res <- llsImpute(object, ...) 

  return(res)
}


##
## This basically copies object@completeObs into the slot
## exprs(ExpressionSet) of an expression set object
##
asExprSet <- function(object, exprSet) {
  if(!inherits(exprSet, "ExpressionSet"))
    stop("Parameter exprSet must be of type ExpressionSet")
  if(!inherits(object, "pcaRes") && !inherits(object, "nniRes"))
    stop("Parameter object must be either of type pcaRes or nniRes")
  if (is.null(completeObs(object)))
    stop("object@completeObs is NULL, exiting")
  if(length(exprs(exprSet)) != length(completeObs(object)))
    stop("Size of exprs(exprSet) and object@completeObs differ. 
Did you really do missing value estimation using this ExpressionSet object?")
  
  exprs(exprSet) <- t(object@completeObs) 
  return(exprSet)
}


plotPcs <- function(object, pcs=1:object@nPcs, type=c("scores", "loadings"), sl=NULL,
                    hotelling=0.95,...) {
  type <- match.arg(type)

  panel <- function(x,y, ...) {
    abline(h=0, v=0, col="black")
    if(!is.null(hotelling)) {
      A <- length(pcs)
      el <- simpleEllipse(x, y, alfa=hotelling)
      lines(el)
    }
    if(is.null(sl))
      points(x, y, ...)
    else
      text(x, y, labels=sl,...)

  }

  switch(type,
         scores = {
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(object@scores[,pcs], labels=labels, panel=panel, upper.panel=NULL,...)
         },
         loadings = {
           if(object@method == "nlpca")
             stop("No loadings plot for Non-linear PCA (the loadings are hidden in a neural network)")
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(object@loadings[,pcs], labels=labels, panel=panel, upper.panel=NULL, ...)
         })
}


setMethod("print", "pcaRes",
          function(x, ...) {
            summary(x)
            cat(x@nVar, "\tVariables\n")
            cat(x@nObs,"\tSamples\n")
            cat(x@missing, "\tNA's\n")
            cat(x@nPcs, "\tCalculated component(s)\n")
            if(x@centered)
              cat("Data was mean centered before running PCA \n")
            else
              cat("Data was NOT mean centered before running PCA \n")
            cat("Scores structure:\n")
            print(dim(x@scores))
            cat("Loadings structure:\n")
            if(x@method == "nlpca") {
              cat("Inverse hierarchical neural network architecture\n")
              cat(drop(x@network@net), "\n")
              cat("Functions in layers\n")
              cat(x@network@fct, "\n")
              cat("hierarchic layer:", x@network@hierarchic$layer, "\n")
              cat("hierarchic coefficients:", x@network@hierarchic$var, "\n")
              cat("scaling factor:", x@network@scalingFactor, "\n")
            }
            else{
              print(dim(x@loadings))
            }
          })

setMethod("leverage", "pcaRes",
          function(object) {
            diag(scores(object) %*%
                 solve(crossprod(scores(object)))  %*% t(scores(object)))
          })

setMethod("show", "pcaRes",
          function(object) {
            summary(object)
            cat(dim(object)["nVar"], "\tVariables\n")
            cat(dim(object)["nObs"],"\tSamples\n")
            cat(object@missing, "\tNA's\n")
            cat(nPcs(object), "\tCalculated component(s)\n")
            if(centered(object))
              cat("Data was mean centered before running PCA \n")
            else
              cat("Data was NOT mean centered before running PCA \n")
            cat("Scores structure:\n")
            print(dim(scores(object)))
            cat("Loadings structure:\n")
            if(method(object) == "nlpca") {
              cat("Inverse hierarchical neural network architecture\n")
              cat(drop(object@network@net), "\n")
              cat("Functions in layers\n")
              cat(object@network@fct, "\n")
              cat("hierarchic layer:", object@network@hierarchic$layer, "\n")
              cat("hierarchic coefficients:", object@network@hierarchic$var, "\n")
              cat("scaling factor:", object@network@scalingFactor, "\n")
            }
            else{
              print(dim(loadings(object)))
            }
          })


## Biplot for pcaRes, uses biplot.default provided by
## package stats. This is basically a copy of the
## biplot.prcomp() function, adapted for a pcaRes object.
## by Kevin Wright
biplot.pcaRes <- function(x, choices=1:2, scale=1, pc.biplot=FALSE, ...) {
  if(length(choices)!=2)
    stop("length of choices must be 2")
  scores <- scores(x)
  n <- nrow(scores)
  lam <- sDev(x)[choices] * sqrt(n)
  if(scale < 0 || scale > 1)
    warning("'scale' is outside [0,1]")
  if(scale != 0) lam <- lam^scale
  else lam <- 1
  if(pc.biplot) lam <- lam/sqrt(n)
  biplot(t(t(scores[,choices])/lam),
         t(t(loadings(x)[, choices]) * lam), , ...)
  invisible()
}
setMethod("biplot", "pcaRes", biplot.pcaRes)

setMethod("print", "nniRes",
          function(x, ...) {
            summary(x)
            cat(dim(x)["nVar"], "\tVariables\n")
            cat(dim(x)["nObs"],"\tSamples\n")
            cat(x@missing, "\tNA's\n")
            cat("k was set to", x@k, "\n")
            if(centered(x))
              cat("Data was mean centered before running LLSimpute \n")
            else
              cat("Data was NOT mean centered before running LLSimpute \n")
          })

setMethod("summary", "pcaRes",
          function(object) {
            if (length(object@R2cum) != length(na.omit(object@R2cum))) {
              warning("Whole row or column was NA in original matrix object. Variance cannot be estimated")
            }
            cat(object@method, "calculated PCA\n")
            cat("Importance of component(s):\n")
            prop <- vector(length=length(object@R2cum), mode="numeric")
            prop[1] <- object@R2cum[1]
            if (length(object@R2cum) > 1) {
              for (i in 2:length(prop)) {
                prop[i] <- object@R2cum[i] - object@R2cum[i-1]
              }
            }
            r <- rbind(prop, object@R2cum)
            rownames(r) <- c("R2", "Cumulative R2")
            colnames(r) <- paste("PC", 1:object@nPcs, sep="")
            print(r, digits=4)
            invisible(r)
          })



predict.pcaRes <- function(object, newdata, pcs=nPcs(object),...) {

  if(!object@method %in% c("ppca", "svd", "nipals", "bpca"))
    stop("predict method not implemented for that type of PCA")

  if(centered(object))  
    newdata <- scale(newdata, object@center, scale=FALSE)
  ## set na's to zero to achieve NIPALS like prediction
  newdata[is.na(newdata)] <- 0
  
  tnew <- newdata %*% loadings(object)[,1:pcs,drop=FALSE]
  xhat <- tcrossprod(tnew,  loadings(object)[,1:pcs,drop=FALSE])
  if(centered(object))  
    xhat <- sweep(xhat, 2, object@center, "+")
  list(scores=tnew, x=xhat)
}

residuals.pcaRes <- function(object, data, nPcs=object@nPcs, ...) {

  data - predict(object, data, nPcs=nPcs)$x
}

fitted.pcaRes <- function(object, data=NULL, nPcs=object@nPcs, ...) {

  ##<..Beg Rdocu..>
  ## ~name~
  ##   fitted.pcaRes
  ## ~title~
  ##   Extract fitted values from PCA.
  ## ~description~
  ##   This function extracts the fitted values from a
  ##   pcaResobject. For PCA methods like SVD, Nipals, PPCA etc this
  ##   is basically just the scores multipled by the loadings,
  ##   forNon-linear PCA the original data is propagated through
  ##   the network to obtain the approximated data.
  ## ~usage~
  ##   fitted.pcaRes(object, data=NULL, nPcs=object@nPcs)
  ## ~arguments~
  ##   ~-object~
  ##     the \code{pcaRes} object of interest.
  ##   ~-data~
  ##     For standard PCA methods this can safely be left null to
  ##     getscores x loadings but if set then the scores are obtained
  ##     byprojecting provided data onto the loadings. 
  ##     If data contains NA values the result will be all NA!!!! (fix this?)
  ##     Non-linear PCA is an exception, here if data is NULL then data is set to
  ##     the completeObs and propaged through the network.
  ##   ~-nPcs~
  ##     The amount of PC's to consider
  ## ~value~
  ##   A matrix with the fitted values.
  ## ~keywords~
  ##   multivariate
  ## ~author~
  ##   Henning Redestig <redestig[at]mpimp-golm.mpg.de>
  ##>..End Rdocu..<

  switch(object@method,
         nlpca = {
           if(is.null(data) & is.null(object@completeObs))
             stop("completeObs slot is empty -- provide the training data")
           if(is.null(data))
             data <- object@completeObs
           data <- t(data)
           if(object@centered)
             data <- sweep(data, 1, object@center)
           recData <- errorHierarchic(object@network, t(object@scores), data)$out[,,nPcs]
           recData <- t(recData / object@network@scalingFactor)
           if(object@centered) 
             recData <- recData + object@center
         },
         {                              #default method
           if(!is.null(data)) {
             if(object@centered)
               data <- sweep(data, 2, object@center)
             scores <- data %*% object@loadings[,1:nPcs, drop=FALSE]
           }
           if(is.null(data))
             scores <- object@scores[,1:nPcs, drop=FALSE]
           recData <- tcrossprod(scores, object@loadings[,1:nPcs, drop=FALSE])
           if(object@centered)
             recData <- t(t(recData) + object@center)
         }
         )
  return(recData)
}

setMethod("fitted", "pcaRes", fitted.pcaRes)


plotR2 <- function(object, nPcs=object@nPcs, type = c("barplot", "lines"),
                   main = deparse(substitute(object)), ...) {
  main <- main    #this is not a typo! the deparse(subsitute(object)) later
                                        ##fails otherwise (dont ask me)
  names(object@sDev) <- paste("PC", 1:nPcs, sep="")
  newx <- list(sdev=object@R2)
  screeplot(newx, nPcs, type, main,...)
}

setMethod("slplot", "pcaRes",
          function(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
                             sl=rownames(object@scores),
                             ll=rownames(object@loadings), hotelling=0.95, rug=TRUE,
                             sub=NULL,...) {

            opar <- par(no.readonly=TRUE)

            cl <- match.call()
            mainArgs <- c(1,match(c("ll", "sl", "scoresLoadings", "sub"),
                                  names(cl), 0))
            scoreArgs <- grep("^s", names(cl)[-mainArgs])
            loadingArgs <- grep("^l", names(cl)[-mainArgs])

            if(!is.null(ll) & length(ll) != nrow(object@loadings))
              stop("Loading labels do not match the object dimensions")
            if(!is.null(sl) & length(sl) != nrow(object@scores))
              stop("Score labels do not match the object dimensions")
            if(is.null(sl))
              sl <- NA
            if(is.null(ll))
              ll <- NA
            
            ## no loadings for non-linear pca
            if(object@method == "nlpca" && scoresLoadings[2])
              scoresLoadings[2] <- FALSE
            
            if(length(pcs) > 2)
              plotPcs(object, pcs, scoresLoadings=scoresLoadings,...)
            else {
              if(is.null(sub))
                sub <- paste(sprintf("%.2f", object@R2cum[min(c(pcs, object@nPcs))]
                                     * 100),
                             "% of the variance explained", sep="")

              if(sum(scoresLoadings) == 2)
                layout(matrix(c(1,2), 1, 2, TRUE), respect=matrix(c(1,1), 1, 2))
              ## exception plot if one dimensional
              if (length(pcs) == 1 || object@nPcs == 1) {
                pcs <- 1
                
                ## score plot
                if(scoresLoadings[1]) {
                  newCall <- call("barplot",
                                  height=object@scores[,pcs],
                                  main="Scores", las=3, ylab=paste("PC", pcs), sub=sub,
                                  names.arg=sl)
                  tmp <- cl[-mainArgs][scoreArgs]
                  names(tmp) <- gsub("^s", "", names(tmp))
                  for(i in 1:length(tmp)) {
                    newCall[[length(newCall) + 1]] <- tmp[[i]]
                    names(newCall)[length(newCall)] <- names(tmp)[i]
                  }
                  eval(newCall)
                }

                ## loadingplot
                if(scoresLoadings[2]) {
                  newCall <- call("barplot",
                                  height=object@loadings[,pcs],
                                  main="Loadings", las=3, ylab=paste("PC", pcs), 
                                  names.arg=ll)
                  if(length(loadingArgs) > 0) {
                    tmp <- cl[-mainArgs][loadingArgs]
                    names(tmp) <- gsub("^l", "", names(tmp))
                    for(i in 1:length(tmp)) {
                      newCall[[length(newCall) + 1]] <- tmp[[i]]
                      names(newCall)[length(newCall)] <- names(tmp)[i]
                    }
                  }
                  eval(newCall)
                }
                return(invisible(TRUE))
              }
              
              ## the score plot
              if(scoresLoadings[1]) {
                ## setup plot
                plotCall <- call("plot",
                                 x=object@scores[,pcs],
                                 main="Scores", ylab=paste("PC", pcs[2]), 
                                 sub=sub, xlab=paste("PC", pcs[1]))
                if(length(scoreArgs) > 0) {
                  tmp <- cl[-mainArgs][scoreArgs]
                  names(tmp) <- gsub("^s", "", names(tmp))
                  for(i in 1:length(tmp)) {
                    plotCall[[length(plotCall) + 1]] <- tmp[[i]]
                    names(plotCall)[length(plotCall)] <- names(tmp)[i]
                  }
                }
                ## add text
                if (!is.null(sl) & !all(is.na(sl))) {
                  plotCall[[length(plotCall) + 1]] <- "n"
                  names(plotCall)[length(plotCall)] <- "type"
                  textCall <- call("text",
                                   x=object@scores[,pcs], labels=sl)
                  if(length(scoreArgs) > 0) {
                    tmp <- cl[-mainArgs][scoreArgs]
                    names(tmp) <- gsub("^s", "", names(tmp))
                    for(i in 1:length(tmp)) {
                      textCall[[length(textCall) + 1]] <- tmp[[i]]
                      names(textCall)[length(textCall)] <- names(tmp)[i]
                    }
                  }
                }
                eval(plotCall)
                if (!is.null(sl) & !all(is.na(sl)))
                  eval(textCall)
                if(rug)
                  rug(object@scores[,1])
                abline(h=0, v=0)
                if(!is.null(hotelling)) {
                  A <- length(pcs)
                  el <- simpleEllipse(object@scores[,pcs[1]],
                                      object@scores[,pcs[2]], alfa=hotelling)
                  lines(el)
                }
              }

              ## the loading plot
              if(scoresLoadings[2]) {
                ## setup plot
                plotCall <- call("plot",
                                 x=object@loadings[,pcs],
                                 main="Loadings", ylab=paste("PC", pcs[2]), 
                                 xlab=paste("PC", pcs[1]))
                if(length(loadingArgs) > 0) {
                  tmp <- cl[-mainArgs][loadingArgs]
                  names(tmp) <- gsub("^l", "", names(tmp))
                  for(i in 1:length(tmp)) {
                    plotCall[[length(plotCall) + 1]] <- tmp[[i]]
                    names(plotCall)[length(plotCall)] <- names(tmp)[i]
                  }
                }
                ## add text
                if (!is.null(ll) & !all(is.na(ll))) {
                  plotCall[[length(plotCall) + 1]] <- "n"
                  names(plotCall)[length(plotCall)] <- "type"
                  textCall <- call("text",
                                   x=object@loadings[,pcs], labels=ll)
                  if(length(loadingArgs) > 0) {
                    tmp <- cl[-mainArgs][loadingArgs]
                    names(tmp) <- gsub("^l", "", names(tmp))
                    for(i in 1:length(tmp)) {
                      textCall[[length(textCall) + 1]] <- tmp[[i]]
                      names(textCall)[length(textCall)] <- names(tmp)[i]
                    }
                  }
                }
                eval(plotCall)
                if (!is.null(ll) & !all(is.na(ll)))
                  eval(textCall)
                abline(h=0, v=0)
              }
            }
            par(opar)
          })


nipalsPca <- function(Matrix, nPcs=2, center = TRUE, completeObs = TRUE,
                      varLimit=1,
                      maxSteps=5000, 
                      threshold=1e-6, verbose=interactive(), ...) {

  ## Convert the object into a matrix (just in case we got a data frame) and do some
  ## basic checks
  Matrix <- as.matrix(Matrix, rownames.force=TRUE)
  if (!checkData(Matrix, verbose = verbose))
    stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

  if (nPcs > ncol(Matrix))
    stop("more components than matrix columns selected, exiting")

  if (center) {
    object <- scale(Matrix, center = TRUE, scale = FALSE)
    means <- attr(object, "scaled:center")
  } else
  object <- Matrix

  missing <- is.na(Matrix)
  nObs <- nrow(object)
  nVar <- ncol(object)

  ##Find a good? starting column -- better way?
  startingColumn <- 1

  ## sum(c(NA, NA), na.rm=TRUE) is 0, but we want NA
  sum.na <- function(x){ ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))}

  TotalSS <- sum(object*object, na.rm=TRUE)

  ph <- rep(0, nVar)
  R2cum <- NULL
  scores <- NULL
  loadings <- NULL
  anotherPc <- TRUE
  l <- 1
  
  while(anotherPc) {
    count <- 0                 #number of iterations done
    th <- object[,startingColumn]   #first column is starting vector for th
    continue <- TRUE
    if(verbose) cat(paste("Calculating PC", l, ": ", sep=""))
    
    while(continue) {
      count <- count+1
      ph <- rep(0, nVar)

      ##Calculate loadings through LS regression
      ##Note: object*th is column-wise multiplication
      tsize <- sum(th * th, na.rm=TRUE)
      ph <- apply(object * (th / tsize), 2, sum.na)
      ##normalize ph based on the available values.
      psize <- sum(ph*ph, na.rm=TRUE)
      ph <- ph / sqrt(psize)
      
      ##Calculate scores through LS regression
      ##Trick: To get row-wise multiplication, use t(object)*ph, then
      ##be sure to use apply(,2,) and NOT apply(,1,)!
      th.old <- th
      th <- apply(t(object) * ph, 2, sum.na)
      
      ##Round up by calculating if convergence condition is met and
      ##checking if it seems to be an neverending loop.
      if (count > maxSteps) {
        stop("Too many iterations, quitting")
      }
      if (t(na.omit(th.old - th)) %*% (na.omit(th.old - th)) <= threshold) {
        continue = FALSE
      }
      if (verbose)cat("*")
    }
    if (verbose) cat(" Done\n")
    object <- object - (th %*% t(ph))
    scores <- cbind(scores, th)
    loadings <- cbind(loadings, ph)
    
    ##cumulative proportion of variance
    R2cum <- cbind(R2cum, 1 - (sum(object*object,na.rm=TRUE) / TotalSS))
    l <- l + 1
    if((!abs(varLimit - 1) < 1e-4 & R2cum[1,l - 1] >= varLimit)| l > nPcs) {
      anotherPc <- FALSE
      nPcs <- l - 1
    }
  }
  R2 <- R2cum[1]
  if(length(R2cum) > 1) R2 <- c(R2, diff(R2cum[1,]))

  if (completeObs) {
    Ye <- scores %*% t(loadings)
    if (center) {
      Ye <- Ye + rep(means, each=nrow(Ye)) # Addition is column-wise
    }
    cObs <- Matrix
    cObs[missing] <- Ye[missing]
  }

  rownames(scores) <- rownames(object)
  colnames(scores) <- paste("PC", 1:nPcs, sep="")
  rownames(loadings) <- colnames(object)
  colnames(loadings) <- paste("PC", 1:nPcs, sep="")
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- cObs
  r@scores <- scores
  r@loadings <- loadings
  r@R2cum <- c(R2cum)
  r@sDev <- apply(scores, 2, sd)
  r@R2 <- R2
  r@nObs <- nObs
  r@nVar <- nVar
  r@varLimit <- varLimit
  r@nPcs <- nPcs
  r@centered <- center
  r@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
  r@method <- "nipals"
  r@missing <- sum(is.na(Matrix))
  return(r)
}


svdPca <- function(Matrix, nPcs=2, center = TRUE, completeObs = FALSE, varLimit=1, verbose=interactive(), ...) {
  
  ## Do some basic checks
  Matrix <- as.matrix(Matrix, rownames.force=TRUE)
  if (!checkData(Matrix, verbose = verbose))
    stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")
  if (nPcs > ncol(Matrix))
    stop("more components than matrix columns selected, exiting")

  if (sum(is.na(Matrix)) > 0)
    stop("SVD PCA cannot handle missing values. Use Nipals PCA, PPCA, BPCA or SVDimpute!")

  if (center) {
    object <- scale(Matrix, center = TRUE, scale = FALSE)
  } else
  object <- Matrix

  pcs <- prcomp(object, center=FALSE, scale.=FALSE, ...)
  imp <- summary(pcs)$importance
  if(varLimit < 1)
    nPcs <- sum(imp[3,] < varLimit) + 1
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- Matrix
  r@scores <- cbind(pcs$x[,1:nPcs])
  colnames(r@scores) <- paste("PC", 1:nPcs, sep = "")
  rownames(r@scores) <- rownames(Matrix) 
  r@loadings <- cbind(pcs$rotation[,1:nPcs])
  colnames(r@loadings) <- paste("PC", 1:nPcs, sep = "")
  rownames(r@loadings) <- colnames(Matrix) 
  r@R2cum <- imp[3,1:nPcs]
  r@sDev <- pcs$sdev[1:nPcs]
  r@R2 <- imp[2,1:nPcs]
  r@nObs <- nrow(object)
  r@nVar <- ncol(object)
  r@varLimit <- varLimit
  r@centered <- center
  r@center <- attr(scale(Matrix, center = TRUE, scale = FALSE), "scaled:center")
  r@nPcs <- nPcs
  r@method <- "svd"
  r@missing <- 0 
  return(r)
}

prep <- function(object, scale=c("none", "pareto", "vector", "UV"), center=TRUE, ...) {
  scale <- match.arg(scale)

  if(inherits(object, "ExpressionSet"))
    object <- t(exprs(object))

  object <- as.matrix(object)
  
  if(center) {
    object <- scale(object, scale=FALSE, center=TRUE)
    cent <- attr(object, "scaled:center")
  }

  if(scale != "none") {
    switch(scale,
           UV = {
             object <- scale(object, scale=TRUE, center=FALSE)
           },
           vector = {
             vectorNorm <- function(y) {
               y <- y / sqrt(sum(y^2))
               return(y)
             }
             object <- t(apply(object, 1, vectorNorm))
           },
           pareto = {
             object <- apply(object, 2,
                             function(object) {
                               return(object / sqrt(sd(object, na.rm=TRUE)))
                             })
           })
  }

  ## recover the center
  if(center)
    attr(object, "scaled:center") <- cent

  attr(object, "scaled") <- scale

  object
}


simpleEllipse <- function(x, y, alfa=0.95, len=200) {

  ##<..Beg Rdocu..>
  ## ~name~
  ##   simpleEllipse
  ## ~title~
  ##   Hotelling's T^2 Ellipse
  ## ~description~
  ##   As described in 'Introduction to multi and megavariate data
  ##   analysis using PCA and PLS' by Eriksson et al. This produces
  ##   very similar ellipse as compared to the ellipse function from
  ##   the ellipse package except that this function assumes that x
  ##   and y are uncorrelated (which they of course are if they are
  ##   scores or loadings from a PCA). This function is mainly
  ##   included to get rid of the dependance on the ellipse package.
  ## ~usage~
  ##   simpleEllipse(x, y, alfa=0.95, len=200)
  ## ~arguments~
  ##   ~-x~
  ##     First PC
  ##   ~-y~
  ##     Second PC
  ##   ~-alfa~
  ##     Significance level of the circle
  ##   ~-len~
  ##     Amount of points in the circle
  ## ~value~
  ##   A matrix with X and Y coordinates for the circle
  ## ~seealso~
  ##   'ellipse'
  ## ~author~
  ##   Henning Redestig <redestig[at]mpimp-golm.mpg.de>
  ##>..End Rdocu..<

  N <- length(x)
  A <- 2
  mypi <- seq(0, 2 * pi, length=len)
  r1 <- sqrt(var(x) * qf(alfa, 2, N - 2) * (2*(N^2 - 1)/(N * (N - 2))))
  r2 <- sqrt(var(y) * qf(alfa, 2, N - 2) * (2*(N^2 - 1)/(N * (N - 2))))
  cbind(r1 * cos(mypi) + mean(x), r2 * sin(mypi) + mean(y))
}


## #############
## some getters

## S3 for those already have been defined as such
scores.pcaRes <- function(object,...) 
  object@scores

loadings.pcaRes <- function(object,...) 
  object@loadings

dim.pcaRes <- function(x)  {
  res <-  c(x@nObs, x@nVar, x@nPcs)
  names(res) <- c("nObs", "nVar", "nPcs")
  res
}

## S4 for new ones
setMethod("nPcs", "pcaRes", function(object, ...) {
  object@nPcs
  })
setMethod("nObs", "pcaRes", function(object, ...) {
  object@nObs
  })
setMethod("nVar", "pcaRes", function(object, ...) {
  object@nVar
  })
setMethod("centered", "pcaRes", function(object, ...) {
  object@centered
  })
setMethod("completeObs", "pcaRes", function(object, ...) {
  object@completeObs
  })
setMethod("method", "pcaRes", function(object, ...) {
  object@method
  })
setMethod("sDev", "pcaRes", function(object, ...) {
  object@sDev
  })

