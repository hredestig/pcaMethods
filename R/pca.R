##' Vector with current valid PCA methods
##' @title List PCA methods
##' @param which the type of methods to get. E.g. only  get the PCA
##' methods based on the classical model where the fitted data is a
##' direct multiplication of scores and loadings. 
##' @return A character vector with the current methods for doing PCA
##' @export 
##' @author Henning Redestig
listPcaMethods <- function(which=c("all", "linear", "nonlinear")) {
  switch(match.arg(which),
         all={
           return(c("svd", "nipals", "rnipals", "bpca", "ppca",
                    "svdImpute", "robustPca", "nlpca", "llsImpute",
                    "llsImputeAll"))
         },
         linear={
           return(c("svd", "nipals", "rnipals", "bpca", "ppca",
                    "svdImpute", "robustPca"))
         },
         nonlinear={
           return("nlpca")
         })
}

##' Perform PCA on a numeric matrix for visualisation, information
##' extraction and missing value imputation.
##'
##' This method is wrapper function for the following set of pca
##' methods:
##'
##' \describe{\item{svd:}{Uses classical \code{prcomp}. See
##' documentation for \code{\link{svdPca}}.}
##' 
##' \item{nipals:}{An iterative method capable of handling small
##' amounts of missing values. See documentation for
##' \code{\link{nipalsPca}}.}
##' 
##' \item{rnipals:}{Same as nipals but implemented in R.}
##' 
##' \item{bpca:}{An iterative method using a Bayesian model to handle
##' missing values. See documentation for \code{\link{bpca}}.}
##' 
##' \item{ppca:}{An iterative method using a probabilistic model to
##' handle missing values. See documentation for \code{\link{ppca}}.}
##' 
##' \item{svdImpute:}{Uses expectation maximation to perform SVD PCA
##' on incomplete data. See documentation for
##' \code{\link{svdImpute}}.}}
##' 
##' Scaling and centering is part of the PCA model and handled by
##' \code{\link{prep}}.
##' @title Perform principal component analysis
##' @param object Numerical matrix with (or an object coercible to
##' such) with samples in rows and variables as columns. Also takes
##' \code{ExpressionSet} in which case the transposed expression
##' matrix is used. Can also be a data frame in which case all
##' numberic variables are used to fit the PCA.
##' @param method One of the methods reported by
##' \code{listPcaMethods()}. Can be left missing in which case the
##' \code{svd} PCA is chosen for data wihout missing values and
##' \code{nipalsPca} for data with missing values
##' @param nPcs Number of principal components to calculate.
##' @param scale Scaling, see \code{\link{prep}}.
##' @param center Centering, see \code{\link{prep}}.
##' @param completeObs Sets the \code{completeObs} slot on the
##' resulting \code{pcaRes} object containing the original data with
##' but with all NAs replaced with the estimates.
##' @param subset A subset of variables to use for calculating the
##' model. Can be column names or indices.
##' @param cv character naming a the type of cross-validation
##' to be performed. 
##' @param ... Arguments to \code{\link{prep}}, the chosen pca
##' method and \code{\link{Q2}}.
##' @return A \code{pcaRes} object.
##' @references
##' Wold, H. (1966) Estimation of principal components and
##' related models by iterative least squares. In Multivariate
##' Analysis (Ed., P.R. Krishnaiah), Academic Press, NY, 391-420.
##'   
##' Shigeyuki Oba, Masa-aki Sato, Ichiro Takemasa, Morito Monden,
##' Ken-ichi Matsubara and Shin Ishii.  A Bayesian missing value
##' estimation method for gene expression profile
##' data. \emph{Bioinformatics, 19(16):2088-2096, Nov 2003}.
##' 
##' Troyanskaya O. and Cantor M. and Sherlock G. and Brown P. and
##' Hastie T. and Tibshirani R. and Botstein D. and Altman RB.  -
##' Missing value estimation methods for DNA microarrays.
##' \emph{Bioinformatics. 2001 Jun;17(6):520-5}.
##' @seealso \code{\link{prcomp}}, \code{\link{princomp}},
##' \code{\link{nipalsPca}}, \code{\link{svdPca}}
##' @examples
##' data(iris)
##' ##  Usually some kind of scaling is appropriate
##' pcIr <- pca(iris, method="svd", nPcs=2)
##' pcIr <- pca(iris, method="nipals", nPcs=3, cv="q2")
##' ## Get a short summary on the calculated model
##' summary(pcIr)
##' plot(pcIr)
##' ## Scores and loadings plot
##' slplot(pcIr, sl=as.character(iris[,5]))
##' 
##' ## use an expressionset and ggplot
##' data(sample.ExpressionSet)
##' pc <- pca(sample.ExpressionSet)
##' df <- merge(scores(pc), pData(sample.ExpressionSet), by=0)
##' library(ggplot2)
##' ggplot(df, aes(PC1, PC2, shape=sex, color=type)) +
##'   geom_point() +
##'   xlab(paste("PC1", pc@R2[1] * 100, "% of the variance")) +
##'   ylab(paste("PC2", pc@R2[2] * 100, "% of the variance"))
##' @export
##' @keywords multivariate
##' @author Wolfram Stacklies, Henning Redestig
pca <- function(object, method, nPcs=2, 
                scale=c("none", "pareto", "vector", "uv"),
                center=TRUE, completeObs=TRUE, subset=NULL,
                cv=c("none","q2"), ...) {
  if(inherits(object, 'data.frame')) {
    num <- vapply(object, is.numeric, logical(1))
    if(sum(num) < 2)
      stop('no numeric data in supplied data.frame')
    Matrix <- as.matrix(object[,num])
  }
  else if(inherits(object, "ExpressionSet")) {
    Matrix <- t(exprs(object))
  } else
  Matrix <- as.matrix(object, rownames.force=TRUE)

  if(!is.null(subset)) 
    Matrix <- Matrix[,subset]
  
  cv <- match.arg(cv)
  scale <- match.arg(scale)

  if (nPcs > ncol(Matrix)) {
    warning("more components than matrix columns requested")
    nPcs <- min(dim(Matrix))
  }
  if (nPcs > nrow(Matrix)) {
    warning("more components than matrix rows requested")
    nPcs <- min(dim(Matrix))
  }

  if (!checkData(Matrix, verbose=interactive()))
    stop("Invalid data format.",
         "Run checkData(data, verbose=TRUE) for details")

  missing <- is.na(Matrix)

  if(missing(method)) {
    if(any(missing))
      method <- 'nipals'
    else
      method <- 'svd'
  }
  if(any(missing) & method == 'svd') {
    warning('data has missing values using nipals instead of user requested svd')
    method <- 'nipals'
  }
  method <- match.arg(method, choices=listPcaMethods())
  
  prepres <- prep(Matrix, scale=scale, center=center, simple=FALSE, ...)

  switch(method,
         svd={
           res <- svdPca(prepres$data, nPcs=nPcs,...)   
         },
         nipals={
           res <- nipalsPca(prepres$data, nPcs=nPcs, ...) 
         },
         rnipals={
           res <- RnipalsPca(prepres$data,  nPcs=nPcs, ...) 
         },
         bpca={
           res <- bpca(prepres$data, nPcs=nPcs, ...) 
         },
         ppca={
           res <- ppca(prepres$data, nPcs=nPcs, ...) 
         },
         svdImpute={
           res <- svdImpute(prepres$data, nPcs=nPcs, ...) 
         },
         robustPca={
           res <- robustPca(prepres$data,  nPcs=nPcs, ...) 
         },
         nlpca={
           res <- nlpca(prepres$data, nPcs=nPcs, ...)
         })

  nPcs <- ncol(res@scores)
  if(is.null(scores(res)) | is.null(loadings(res)) |
     is.null(R2cum(res)) | is.null(method(res)))
    stop(paste("bad result from pca method", method))

  colnames(res@scores) <- paste("PC", 1:nPcs, sep="")
  rownames(res@scores) <- rownames(Matrix)
  if(all(dim(loadings(res)) == c(ncol(Matrix), nPcs))) {
    colnames(res@loadings) <- paste("PC", 1:nPcs, sep="")
    rownames(res@loadings) <- colnames(Matrix)
  }
  if(!is.null(subset))
    res@subset <- subset
  res@missing <- missing
  res@nPcs <- nPcs
  res@nObs <- nrow(Matrix)
  res@nVar <- ncol(Matrix)
  res@sDev <- apply(scores(res), 2, sd)
  res@center <- prepres$center
  res@centered <- center
  res@scale <- prepres$scale
  res@scaled <- scale 
  res@R2 <- res@R2cum[1]
  if(length(res@R2cum) > 1)
    res@R2 <- c(res@R2, diff(res@R2cum))

  if (completeObs) {
    cObs <- Matrix
    if(method %in% listPcaMethods("nonlinear"))
      cObs[missing] <- fitted(res, Matrix, pre=TRUE, post=TRUE)[missing]
    else
      cObs[missing] <- fitted(res, post=TRUE)[missing]
    res@completeObs <- cObs
  }
  if(cv == "q2")
    res@cvstat <- Q2(res, Matrix, nruncv=1, ...)

  return(res)
}

##' Wrapper function for imputation methods based on nearest neighbour
##' clustering. Currently llsImpute only.
##'
##' This method is wrapper function to llsImpute, See documentation
##' for \code{link{llsImpute}}.
##' @title Nearest neighbour imputation
##' @param object Numerical matrix with (or an object coercible to
##' such) with samples in rows and variables as columns. Also takes
##' \code{ExpressionSet} in which case the transposed expression
##' matrix is used.
##' @param method For convenience one can pass a large matrix but only
##' use the variable specified as subset. Can be colnames or indices.
##' @param subset Currently "llsImpute" only.
##' @param ... Further arguments to the chosen method.
##' @return A \code{clusterRes} object. Or a list containing a
##' clusterRes object as first and an ExpressionSet object as second
##' entry if the input was of type ExpressionSet.
##' @export
##' @seealso \code{\link{llsImpute}}, \code{\link{pca}}
##' @keywords multivariate
##' @examples
##' data(metaboliteData)
##' llsRes <- nni(metaboliteData, k=6, method="llsImpute", allGenes=TRUE)
##' @author Wolfram Stacklies
nni <- function(object, method=c("llsImpute"), subset=numeric(), ...) {
  isExprSet <- FALSE
  if(inherits(object, "ExpressionSet")) {
    set <- object
    isExprSet <- TRUE
    object <- t(exprs(object))
  }

  method <- match.arg(method)

  if ( !checkData(as.matrix(object), verbose=interactive()) )
    stop("Invalid data format, exiting...\n",
         "Run checkData(data, verbose=TRUE) for details\n")

  missing <- sum(is.na(object))
  if(length(subset) > 0)
    object <- object[,subset]

  res <- llsImpute(object, ...) 

  return(res)
}

##' A function that can be used to visualise many PCs plotted against
##' each other
##'
##' Uses \code{\link{pairs}} to provide side-by-side plots. Note that
##' this function only plots scores or loadings but not both in the
##' same plot.
##' @title Plot many side by side scores XOR loadings plots
##' @param object \code{pcaRes} a pcaRes object
##' @param pcs \code{numeric} which pcs to plot
##' @param type \code{character} Either "scores" or "loadings" for
##' scores or loadings plot respectively
##' @param sl \code{character} Text labels to plot instead of a point,
##' if NULL points are plotted instead of text
##' @param hotelling \code{numeric} Significance level for the
##' confidence ellipse. NULL means that no ellipse is drawn.
##' @param ... Further arguments to \code{\link{pairs}} on which this
##' function is based. 
##' @return None, used for side effect.
##' @seealso \code{prcomp}, \code{pca}, \code{princomp}, \code{slplot}
##' @export
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4], nPcs=3,  method="svd")
##' plotPcs(pcIr, col=as.integer(iris[,4]) + 1)
##' @keywords multivariate
##' @author Henning Redestig
plotPcs <- function(object,
                    pcs=1:nP(object), type=c("scores", "loadings"), sl=NULL,
                    hotelling=0.95,  ...) {
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
         scores={
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(scores(object)[,pcs], labels=labels,
                 panel=panel, upper.panel=NULL,...)
         },
         loadings={
           if(method(object) == "nlpca")
             stop("Loadings plot not applicable for non-linear PCA")
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(loadings(object)[,pcs], labels=labels, panel=panel,
                 upper.panel=NULL, ...)
         })
}


##' A wrapper function for \code{prcomp} to deliver the result as a
##' \code{pcaRes} method. Supplied for  compatibility with  the rest
##' of the pcaMethods package.  It is not recommended to use this
##' function directely but rather to use  the \code{pca()} wrapper
##' function. 
##' @title Perform principal component analysis using singular value
##' decomposition
##' @param Matrix Pre-processed (centered and possibly scaled)
##' numerical matrix samples in rows and variables as columns. No
##' missing values allowed.
##' @param nPcs Number of components that should be extracted.
##' @param varLimit Optionally the ratio of variance that should be
##' explained. \code{nPcs} is ignored if varLimit < 1
##' @param verbose Verbose complaints to matrix structure
##' @param ... Only used for passing through arguments.
##' @return A \code{pcaRes} object.
##' @seealso \code{prcomp}, \code{princomp}, \code{pca}
##' @examples
##' data(metaboliteDataComplete)
##' mat <- prep(t(metaboliteDataComplete))
##' pc <- svdPca(mat, nPcs=2)
##' ## better use pca()
##' pc <- pca(t(metaboliteDataComplete), method="svd", nPcs=2)
##' \dontshow{stopifnot(sum((fitted(pc) - t(metaboliteDataComplete))^2, na.rm=TRUE) < 200)}
##' @export
##' @keywords multivariate
##' @author Henning Redestig
svdPca <- function(Matrix, nPcs=2, 
                   varLimit=1, verbose=interactive(), ...) {

  pcs <- prcomp(Matrix, center=FALSE, scale.=FALSE)
  imp <- summary(pcs)$importance
  if(varLimit < 1)
    nPcs <- sum(imp[3,] < varLimit) + 1

  res <- new("pcaRes")
  res@scores <- cbind(pcs$x[,1:nPcs])
  res@loadings <- cbind(pcs$rotation[,1:nPcs])
  res@R2cum <- imp[3,1:nPcs]
  res@varLimit <- varLimit
  res@method <- "svd"
  return(res)
}

##' Get a confidence ellipse for uncorrelated bivariate data
##'
##' As described in 'Introduction to multi and megavariate data analysis
##' using PCA and
##' PLS' by Eriksson et al. This produces very similar ellipse as
##' compared to the ellipse function the ellipse package except that
##' this function assumes that and y are uncorrelated (which they of
##' are if they are scores or loadings from a PCA). 
##' @title Hotelling's T^2 Ellipse
##' @param x first variable
##' @param y second variable
##' @param alfa confidence level of the circle
##' @param len Number of points in the circle
##' @seealso ellipse
##' @author Henning Redestig
##' @return A matrix with X and Y coordinates for the circle
simpleEllipse <- function(x, y, alfa=0.95, len=200) {
  N <- length(x)
  A <- 2
  mypi <- seq(0, 2 * pi, length=len)
  r1 <- sqrt(var(x) * qf(alfa, 2, N - 2) * (2*(N^2 - 1)/(N * (N - 2))))
  r2 <- sqrt(var(y) * qf(alfa, 2, N - 2) * (2*(N^2 - 1)/(N * (N - 2))))
  cbind(r1 * cos(mypi) + mean(x), r2 * sin(mypi) + mean(y))
}


# .onLoad <- function(libname, pkgname) {
#   require("methods")
# }
