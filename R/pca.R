pca <- function(object, method=c("svd", "nipals", "bpca", "ppca", "svdImpute"), subset=numeric(),...) {
  
  if(inherits(object, "exprSet")) 
    object <- t(exprs(object))

  method <- match.arg(method)

  # Do some basic checks of the data. We exit if the data contains NaN or Inf
  # values or is not numeric.
  if ( !checkData(as.matrix(object), verbose = interactive()) )
	stop("Invalid data format, exiting...\n",
	     "Run checkData(data, verbose = TRUE) for details\n")
	
  missing <- sum(is.na(object))

  if(length(subset) > 0)
    object <- object[,subset]

  #object <- prep(object, ...)

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
         })

  scaled = attr(object, "scaled")
  if (!is.null(scaled)) {
    res@scaled <- scaled
  } else
    res@scaled <- "none"

  if(length(subset) > 0)
    res@subset <- subset
    
  return(res)
}

plotPcs <- function(object, pc=1:object@nPcs, scoresLoadings=c(TRUE, FALSE),...) {

  ## number of plot areas needed
  ba <- (length(pc)*(length(pc) - 1)) / 2
  mf <- c(floor(sqrt(ba)), ceiling(sqrt(ba)))
  if(mf[1] * mf[2] < ba)
    mf[1] <- mf[1] + 1
  mat <- matrix(c(1:ba, rep(0, mf[1]*mf[2] - ba)), ncol=mf[1], nrow=mf[2], byrow=TRUE)
  layout(mat)

  perm <- function(object) {
    res <- NULL
    tri <- matrix(TRUE, length(object), length(object))
    tri[upper.tri(tri, TRUE)] <- FALSE
    for(i in 1:length(object))
      for(j in 1:length(object))
        if(tri[i,j])
          res <- rbind(res, c(object[i],object[j]))
    res
  }

  pp <- t(apply(perm(pc), 1, sort))

  for(i in 1:nrow(pp))
    slplot(object, pcs=pp[i,], scoresLoadings=scoresLoadings, rug=FALSE, sub="",...)
}


setMethod("print", "pcaRes",
          function(x, ...) {
            summary(x)
            cat(x@nVar, "\tVariables\n")
            cat(x@nObs,"\tSamples\n")
            cat(x@missing, "\tNA's\n")
            cat(x@nPcs, "\tCalculated component(s)\n")
            cat("*** Data was ")
            if(x@scaled != "none") { cat("scaled (using ",x@scaled,") ", sep="") }
            else { cat("NOT scaled ") }
            if(x@centered && (x@scaled == "none")) { cat("but centered ***\n") }
            else if(x@centered && (x@scaled != "none")) { cat("and centered ***\n") }
            else if(!x@centered && (x@scaled != "none")) { cat("but NOT centered ***\n") }
            else if(!x@centered && (x@scaled == "none")) { cat("and NOT centered ***\n") }
            cat("Scores structure:\n")
            print(dim(x@scores))
            cat("Loadings structure:\n")
            print(dim(x@loadings))
          })

setMethod("summary", "pcaRes",
          function(object) {
            if (length(object@R2cum) != length(na.omit(object@R2cum))) {
              warning("Whole row or column was NA in original matriobject. Variance cannot be estimated")
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
            rownames(r) <- c("Proportion of Variance", "Cumulative Proportion")
            colnames(r) <- paste("PC", 1:object@nPcs, sep="")
            print(r, digits=4)
            invisible(r)
          })

setMethod("screeplot", "pcaRes",
          function(x, npcs=min(10, length(x@R2)), type = c("barplot", "lines"),
                   main = deparse(substitute(x)), ...) {
            newx <- list(sDev=x@sDev)
            screeplot(newx, nPcs, type, main,...)
          })

setMethod("slplot", "pcaRes",
          function(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE), sl=rownames(object@scores),
                   ll=rownames(object@loadings), hotelling=0.95, rug=TRUE,sub=NULL,...) {
            if(length(pcs) > 2)
              plotPcs(object, pcs, scoresLoadings=scoresLoadings,...)

            else {
              if(is.null(sub))
                sub <- paste(sprintf("%.2f", object@R2cum[pcs] * 100),
                             "% of the variance explained", sep="")

              if(sum(scoresLoadings) == 2)
                layout(matrix(c(1,2), 1, 2, TRUE), respect=matrix(c(1,1), 1, 2))
              if (length(pcs) == 1 || object@nPcs == 1) {
                pcs <- 1
                if(scoresLoadings[1]) {
                  barplot(object@scores[,pcs], main="Scores",
                          sub=sub, las=3, names.arg=sl,
                          ylab=paste("PC", pcs), ...)
                }
                if(scoresLoadings[2])
                  barplot(object@loadings[,pcs], main="Loadings", ylab=paste("PC", pcs), las=3, ...)
                return(TRUE)
              }
              if(scoresLoadings[1]) {
                if (!is.null(sl)) {
                  plot(object@scores[,pcs], type="n", main="Scores",
                       sub=sub, ylab=paste("PC", pcs[2]),
                       xlab=paste("PC", pcs[1]),...)
                  text(object@scores[,pcs], sl,...)
                  if(rug)
                    rug(object@scores[,1])
                }
                else {
                  plot(object@scores[,pcs], main="Scores",
                       sub=sub,ylab=paste("PC", pcs[2]),
                       xlab=paste("PC", pcs[1]),...)
                }
                abline(h=0, v=0)
                if(!is.null(hotelling)) {
                  A <- length(pcs)
                  lines(ellipse(cov(object@scores[,pcs]), level=hotelling))
                }
              }
              if(scoresLoadings[2]) {
                if (!is.null(ll)) {
                  plot(object@loadings[,pcs], type="n", main="Loadings", ylab=paste("PC", pcs[2]),
                       xlab=paste("PC", pcs[1]), ...)
                  text(object@loadings[,pcs], ll,...)
                }
                else {
                  plot(object@loadings[,pcs], main="Loadings",ylab=paste("PC", pcs[1]), ...)
                }
                abline(h=0, v=0)
              }
            }
          })

  
nipalsPca <- function(Matrix, nPcs=2, center = TRUE, completeObs = TRUE, varLimit=1, maxSteps=5000, 
                      threshold=1e-6, verbose=interactive(), ...) {

  ## Convert the object into a matrix (just in case we got a data frame) and do some
  ## basic checks
  Matrix <- as.matrix(Matrix)
  if (!checkData(Matrix, verbose = verbose))
    stop("Invalid data format! Use checkData(Matrix, verbose = TRUE) for details.\n")

  if (nPcs > ncol(Matrix))
    stop("more components than matrix columns selected, exiting")

  if (center) {
    object <- scale(Matrix, center = TRUE)
    means <- attr(object, "scaled:center")
  } else
    object <- Matrix

  missing <- is.na(Matrix)
  nObs <- nrow(object)
  nVar <- ncol(object)

  ##Find a good starting column
  ##(one where not all values are NA) and do some initial checking
  for (i in 1:nVar) {         
    if (length(na.omit(object[,i])) != 0) {
      startingColumn <- i
    }
    else {
      warning("Whole column was found NA in original matrix")
    }
  }
  for (i in 1:nObs) {
    if (length(na.omit(object[i,])) == 0) {
      warning("Whole row was found NA in original matrix")
    }
  }

  scaledobject <- object

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
      for (i in 1:nVar) {
        ph[i] <- sum(object[,i] * th, na.rm=TRUE) / sum(th^2, na.rm=TRUE)
        ##if the whole column is NA the loading of that variable should also be NA
        if (length(na.omit(object[,i])) == 0) { 
          ph[i] <- NA
        }
      }
      ##normalize ph based on the available values.
      ph <- ph / sqrt(sum(na.omit(ph)^2))

      ##Calculate scores through LS regression
      th.old <- th
      th <- rep(0, nObs)
      for (i in 1:nObs) {
        th[i] <- sum(object[i,] * ph, na.rm=TRUE) / sum(ph^2, na.rm=TRUE)
        if (length(na.omit(object[i,])) == 0) {
          th[i] <- NA
        }
      }
      
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
    R2cum <- cbind(R2cum, 1 - (sum(object^2,na.rm=TRUE) / sum(scaledobject^2, na.rm=TRUE)))
    l <- l + 1
    if (R2cum[1,l - 1] >= varLimit || l > nPcs) {
      anotherPc <- FALSE
      nPcs <- l - 1
    }
  }
  R2 <- vector(length=length(R2cum), mode="numeric")
  R2[1] <- R2cum[1]
  if (length(R2cum) > 1) {
    for (i in 2:length(R2)) {
      R2[i] <- R2cum[i] - R2cum[i-1]
    }
  }

  if (completeObs) {
    Ye <- scores %*% t(loadings)
    if (center) {
      for(i in 1:ncol(Ye)) {
        Ye[,i] <- Ye[,i] + means[i]
      }
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
  r@scaled <- "none"
  r@centered <- center
  r@center <- attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center")
  r@method <- "nipals"
  r@missing <- sum(is.na(Matrix))
  return(r)
}


svdPca <- function(Matrix, nPcs=2, center = TRUE, completeObs = FALSE, varLimit=1, ...) {
  
  ## Do some basic checks
  Matrix <- as.matrix(Matrix)
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

  pcs <- prcomp(object, center=FALSE, scale.=FALSE)
  imp <- summary(pcs)$importance
  if(varLimit < 1)
    nPcs <- sum(imp[3,] < varLimit) + 1
  r <- new("pcaRes")
  if (completeObs)
    r@completeObs <- Matrix
  r@center <- attr(scale(Matrix, center=TRUE, scale=FALSE), "scaled:center")
  r@scores <- cbind(pcs$x[,1:nPcs])
  r@loadings <- cbind(pcs$rotation[,1:nPcs])
  r@R2cum <- imp[3,1:nPcs]
  r@sDev <- pcs$sdev[1:nPcs]
  r@R2 <- imp[2,1:nPcs]
  r@nObs <- nrow(object)
  r@nVar <- ncol(object)
  r@varLimit <- varLimit
  r@scaled <- "none"
  r@centered <- center
  r@nPcs <- nPcs
  r@method <- "svd"
  r@missing <- 0 
  return(r)
}

prep <- function(object, scale=c("none", "pareto", "vector", "UV"), center=TRUE, ...) {
  scale <- match.arg(scale)

  if(inherits(object, "exprSet"))
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

