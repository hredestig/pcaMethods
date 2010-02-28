## Hope the fixes a strange bug during package test on bioconductor
.onLoad <- function(libname, pkgname) {
  require("methods")
}

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
         all = {
           return(c("svd", "nipals", "rnipals", "bpca", "ppca",
                    "svdImpute", "robustPca", "nlpca"))
         },
         linear = {
           return(c("svd", "nipals", "rnipals", "bpca", "ppca",
                    "svdImpute", "robustPca"))
         },
         nonlinear = {
           return("nlpca")
         })
}

##' Can be used for computing PCA on a numeric matrix for
##' visualisation, information extraction and missing value
##' imputation.
##'
##' This method is wrapper function for the following set of pca
##' methods:
##'
##' \describe{\item{svd:}{Uses classical \code{prcomp}.See
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
##' matrix is used.
##' @param method One of the methods reported by \code{pcaMethods()}
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
##' pcIr <- pca(iris[,1:4], method="svd", nPcs=2)
##' pcIr <- pca(iris[,1:4], method="nipals", nPcs=3, cv="q2")
##' ## Get a short summary on the calculated model
##' summary(pcIr)
##' plot(pcIr)
##' ## Scores and loadings plot
##' slplot(pcIr, sl=as.character(iris[,5]))
##' @export
##' @keywords multivariate
##' @author Wolfram Stacklies, Henning Redestig
pca <- function(object, method=listPcaMethods(), nPcs=2, 
                scale=c("none", "pareto", "vector", "uv"),
                center=TRUE, completeObs=TRUE, subset=NULL,
                cv=c("none","q2"), ...) {
  
  if(inherits(object, "ExpressionSet")) {
    Matrix <- t(exprs(object))
  } else
  Matrix <- as.matrix(object, rownames.force=TRUE)

  method <- match.arg(method)
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

  if(!is.null(subset)) 
    Matrix <- Matrix[,subset]
  
  missing <- is.na(Matrix)

  if(any(missing) & method == "svd") 
    method <- "nipals"
  
  prepres <- prep(Matrix, scale=scale, center=center, simple=FALSE, ...)

  switch(method,
         svd = {
           res <- svdPca(prepres$data, nPcs=nPcs,...)   
         },
         nipals = {
           res <- nipalsPca(prepres$data, nPcs=nPcs, ...) 
         },
         rnipals = {
           res <- RnipalsPca(prepres$data,  nPcs=nPcs, ...) 
         },
         bpca = {
           res <- bpca(prepres$data, nPcs=nPcs, ...) 
         },
         ppca = {
           res <- ppca(prepres$data, nPcs=nPcs, ...) 
         },
         svdImpute = {
           res <- svdImpute(prepres$data, nPcs=nPcs, ...) 
         },
         robustPca = {
           res <- robustPca(prepres$data,  nPcs=nPcs, ...) 
         },
         nlpca = {
           res <- nlpca(prepres$data, nPcs=nPcs, ...)
         })

  if(is.null(scores(res)) | is.null(loadings(res)) |
     is.null(R2cum(res)) | is.null(method(res)))
    stop(paste("bad result from pca method", method))

  colnames(res@scores) <- paste("PC", 1:nPcs, sep = "")
  rownames(res@scores) <- rownames(Matrix)
  if(all(dim(loadings(res)) == c(ncol(Matrix), nPcs))) {
    colnames(res@loadings) <- paste("PC", 1:nPcs, sep = "")
    rownames(res@loadings) <- colnames(Matrix)
  }
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

  if ( !checkData(as.matrix(object), verbose = interactive()) )
    stop("Invalid data format, exiting...\n",
         "Run checkData(data, verbose = TRUE) for details\n")

  missing <- sum(is.na(object))
  if(length(subset) > 0)
    object <- object[,subset]

  res <- llsImpute(object, ...) 

  return(res)
}

##' This function can be used to conveniently replace the expression
##' matrix in an \code{ExpressionSet} with the completed data from a
##' \code{pcaRes} object.
##'
##' This is not a standard \code{as} function as \code{pcaRes}
##' object alone not can be converted to an \code{ExpressionSet} (the
##' \code{pcaRes} object does not hold any \code{phenoData} for
##' example).
##' @title Convert pcaRes object to an expression set
##' @param object \code{pcaRes} -- The object containing the completed
##' data.
##' @param exprSet \code{ExpressionSet} -- The object passed on to
##' \code{pca} for missing value estimation.  
##' @return An object without missing values of class \code{ExpressionSet}.
##' @export 
##' @author Wolfram Stacklies \cr CAS-MPG Partner Institute for
##' Computational Biology, Shanghai, China
##' @keywords multivariate
asExprSet <- function(object, exprSet) {
  if(!inherits(exprSet, "ExpressionSet"))
    stop("Parameter exprSet must be of type ExpressionSet")
  if(!inherits(object, "pcaRes") & !inherits(object, "nniRes"))
    stop("Parameter object must be either of type pcaRes or nniRes")
  if (is.null(completeObs(object)))
    stop("completeObs(object) is NULL, exiting")
  if(!all(dim(exprs(exprSet)) == dim(t(completeObs(object)))))
    stop("Dimensions of exprs(exprSet) and completeObs(object) do not match. 
Did you really do missing value estimation using this ExpressionSet object?")
  exprs(exprSet) <- t(completeObs(object))
  return(exprSet)
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
         scores = {
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(scores(object)[,pcs], labels=labels,
                 panel=panel, upper.panel=NULL,...)
         },
         loadings = {
           if(method(object) == "nlpca")
             stop("Loadings plot not applicable for non-linear PCA")
           labels <- paste("PC", pcs, "\n", "R^2 =", round(object@R2[pcs], 2))
           pairs(loadings(object)[,pcs], labels=labels, panel=panel,
                 upper.panel=NULL, ...)
         })
}

##' Print basic information about pcaRes object
##' @title Print/Show for pcaRes
##' @param x a pcaRes object
##' @param ...  not used
##' @return nothing, used for its side effect
##' @export
##' @author Henning Redestig
showPcaRes <- function(x, ...) {
  summary(x)
  cat(nVar(x), "\tVariables\n")
  cat(nObs(x),"\tSamples\n")
  cat(nmissing(x), "\tNAs (",
      round(100 * nmissing(x) / (nObs(x) * nVar(x)),
            getOption("str")$digits.d), "%)\n")
  cat(nP(x), "\tCalculated component(s)\n")
  if(centered(x))
    cat("Data was mean centered before running PCA \n")
  else
    cat("Data was NOT mean centered before running PCA \n")
  if(scaled(x))
    cat("Data was scaled before running PCA \n")
  else
    cat("Data was NOT scaled before running PCA \n")
  cat("Scores structure:\n")
  print(dim(scores(x)))
  cat("Loadings structure:\n")
  if(method(x) == "nlpca") {
    cat("Inverse hierarchical neural network architecture\n")
    cat(drop(x@network@net), "\n")
    cat("Functions in layers\n")
    cat(x@network@fct, "\n")
    cat("hierarchic layer:", x@network@hierarchic$layer, "\n")
    cat("hierarchic coefficients:", x@network@hierarchic$var, "\n")
    cat("scaling factor:", x@network@scalingFactor, "\n")
  }
  else{
    print(dim(loadings(x)))
  }
}

##' Print a brief description of nniRes model
##' @title Print a nniRes model
##' @param x An \code{nniRes} object
##' @param ... Not used
##' @return Nothing, used for side-effect
##' @export
##' @author Henning Redestig
showNniRes <- function(x, ...) {
  summary(x)
  cat(dim(x)["nVar"], "\tVariables\n")
  cat(dim(x)["nObs"],"\tSamples\n")
  cat(nmissing(x), "\tNAs (",
      round(100 * nmissing(x) / (nObs(x) * nVar(x)),
            getOption("str")$digits.d), "%)\n")
  cat("k was set to", x@k, "\n")
  if(centered(x))
    cat("Data was mean centered before running LLSimpute \n")
  else
    cat("Data was NOT mean centered before running LLSimpute \n")
  if(scaled(x))
    cat("Data was scaled before running LLSimpute \n")
  else
    cat("Data was NOT scaled before running LLSimpute \n")
}
setMethod("print", "nniRes", showNniRes)
##' Print basic info
##' @seealso \code{\link{showPcaRes}} \code{\link{showNniRes}}
##' @exportMethod print
##' @aliases print,pcaRes-method print,nniRes-method
setMethod("print", "pcaRes", showPcaRes)
setMethod("show", "nniRes", function(object) showNniRes(object))
##' Show pcaRes / nniRes objects.
##' @seealso \code{\link{showPcaRes}} \code{\link{showNniRes}}
##' @importFrom methods,show
##' @exportMethod show
##' @aliases show,pcaRes-method show,nniRes-method
setMethod("show", "pcaRes", function(object) showPcaRes(object))

##' Visualize two-components simultaneously
##'
##' This is a method for the generic function 'biplot'.  There is
##' considerable confusion over the precise definitions: those of the
##' original paper, Gabriel (1971), are followed here.  Gabriel and
##' Odoroff (1990) use the same definitions, but their plots actually
##' correspond to \code{pc.biplot = TRUE}. 
##' @title Plot a overlaid scores and loadings plot
##' @param x a pcaRes object
##' @param choices which two pcs to plot 
##' @param scale The variables are scaled by
##' \eqn{\lambda^{scale}}{lambda^scale} and the observations are
##' scaled    by \eqn{\lambda^{scale}}{lambda ^ (1-scale)} where
##' \code{lambda} are  the singular values as computed by
##' \code{princomp}.  Normally  \eqn{0\le{}scale\le{}1}{0 <= scale <=
##' 1}, and a warning will be issued if the specified 'scale' is
##' outside this range.
##' @param pc.biplot If true, use what Gabriel (1971) refers to as a
##' "principal component biplot", with \eqn{\lambda=1}{lambda = 1} and
##' observations scaled up by sqrt(n) and variables scaled down by
##' sqrt(n). Then the inner products between variables approximate
##' covariances and distances between observations approximate
##' Mahalanobis distance. 
##' @param ... optional arguments to be passed to
##' \code{biplot.default}. 
##' @return a plot is produced on the current graphics device.
##' @aliases biplot.pcaRes
##' @export
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4])
##' biplot(pcIr)
##' @seealso \code{prcomp}, \code{pca}, \code{princomp}
##' @author Kevin Wright, Adapted from \code{biplot.prcomp}
##' @keywords multivariate
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
##' Biplot for pcaRes method. 
##' @seealso \code{\link{biplot.pcaRes}}
##' @aliases biplot,pcaRes-method
##' @importFrom stats,biplot
setMethod("biplot", "pcaRes", biplot.pcaRes)

##' Print a brief description of the PCA model
##' @title Summary of PCA model
##' @param object a pcaRes object
##' @param ... Not available
##' @usage summary(object, ...)
##' @return Nothing, used for side-effect
##' @exportMethod summary
##' @aliases summary summary,pcaRes-method
##' @author Henning Redestig
setMethod("summary", "pcaRes",
          function(object) {
            cat(method(object), "calculated PCA\n")
            cat("Importance of component(s):\n")
            prop <- vector(length=length(R2cum(object)), mode="numeric")
            prop[1] <- R2cum(object)[1]
            if (length(R2cum(object)) > 1) {
              for (i in 2:length(prop)) {
                prop[i] <- R2cum(object)[i] - R2cum(object)[i-1]
              }
            }
            r <- rbind(prop, R2cum(object))
            rownames(r) <- c("R2", "Cumulative R2")
            colnames(r) <- paste("PC", 1:nP(object), sep="")
            print(r, digits=4)
            invisible(r)
          })

##' Predict data using PCA model
##'
##' This function extracts the predict values from a pcaRes object for
##' the PCA methods SVD, Nipals, PPCA and BPCA.  Newdata is first
##' centered if the PCA model was and then scores (\eqn{T}) and data
##' (\eqn{X}) is 'predicted' according to :
##' \eqn{\hat{T}=X_{new}P}{That=XnewP}
##' \eqn{\hat{X}_{new}=\hat{T}P'}{Xhat=ThatP'}.  Missing values are
##' set to zero before matrix multiplication to achieve NIPALS like
##' treatment of missing values. 
##' @title Predict values from PCA.
##' @param object \code{pcaRes} the \code{pcaRes} object of interest.
##' @param newdata \code{matrix} new data with same number of columns
##' as the used to compute \code{object}.
##' @param pcs \code{numeric} The number of PC's to consider
##' @param pre pre-process \code{newdata} based on the pre-processing
##' chosen for the PCA model
##' @param post unpre-process the final data (add the center back etc)
##' @param ... Not passed on anywhere, included for S3 consistency.
##' @return A list with the following components: \item{scores}{The
##' predicted scores} \item{x}{The predicted data}
##' @aliases predict.pcaRes 
##' @keywords multivariate
##' @examples
##' data(iris)
##' hidden <- sample(nrow(iris), 50)
##' pcIr <- pca(iris[-hidden,1:4])
##' pcFull <- pca(iris[,1:4])
##' irisHat <- predict(pcIr, iris[hidden,1:4])
##' cor(irisHat$scores[,1], scores(pcFull)[hidden,1])
##' @export 
##' @author Henning Redestig
predict.pcaRes <- function(object, newdata, pcs=nP(object),
                           pre=TRUE, post=TRUE, ...) {
  if(!method(object) %in% listPcaMethods("linear"))
    stop("predict method not implemented for that type of PCA")
  if(pre)
    newdata <- prep(newdata, scl(object), center(object))
  ## set na's to zero to achieve NIPALS like prediction
  newdata[is.na(newdata)] <- 0
  tnew <- newdata %*% loadings(object)[,1:pcs,drop=FALSE]
  xhat <- tcrossprod(tnew,  loadings(object)[,1:pcs,drop=FALSE])
  if(post)
    xhat <- prep(xhat, scl(object), center(object), reverse=TRUE)
  list(scores=tnew, x=xhat)
}
##' Predict PCA data.
##' @seealso \code{\link{predict.pcaRes}}
##' @importFrom stats,predict
##' @exportMethod predict
##' @aliases predict,pcaRes-method
setMethod("predict", "pcaRes", predict.pcaRes)

##' This function extracts the residuals values from a pcaRes object
##' for the PCA methods SVD, Nipals, PPCA and BPCA
##' @title Residuals values from a PCA model.
##' @param object \code{pcaRes} the \code{pcaRes} object of interest.
##' @param data \code{matrix} The data that was used to calculate the
##' PCA model (or a different dataset to e.g. adress its proximity to
##' the model). 
##' @param ... Passed on to \code{\link{predict.pcaRes}}. E.g. setting
##' the number of used components.
##' @return A \code{matrix} with the residuals
##' @aliases residuals.pcaRes
##' @keywords multivariate
##' @export
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4])
##' head(residuals(pcIr, iris[,1:4]))
##' @author Henning Redestig
residuals.pcaRes <- function(object, data=completeObs(object), ...) {
  if(is.null(data))
    stop("data missing when calculating residuals")
  data - predict(object, data, ...)$x
}
##' Residuals of PCA data.
##' @seealso \code{\link{residuals.pcaRes}}
##' @importFrom stats,residuals
##' @exportMethod residuals
##' @aliases residuals,pcaRes-method
setMethod("residuals", "pcaRes", residuals.pcaRes)
##' Residuals of PCA data.
##' @seealso \code{\link{residuals.pcaRes}}
##' @importFrom stats,resid
##' @exportMethod resid
##' @aliases resid,pcaRes-method
setMethod("resid", "pcaRes", residuals.pcaRes)

##' Fitted values of a PCA model
##' 
##' This function extracts the fitted values from a pcaResobject. For
##' PCA methods like SVD, Nipals, PPCA etc this is basically just the
##' scores multipled by the loadings and adjusted for pre-processing.
##' for non-linear PCA the original data is propagated through the
##' network to obtain the approximated data.
##' @title Extract fitted values from PCA.
##' @param object the \code{pcaRes} object of interest.
##' @param data For standard PCA methods this can safely be left null
##' to get scores x loadings but if set, then the scores are obtained
##' by projecting provided data onto the loadings.  If data contains
##' missing values the result will be all NA. Non-linear PCA is an
##' exception, here if data is NULL then data is set to the
##' completeObs and propaged through the network.
##' @param nPcs The number of PC's to consider
##' @param pre pre-process \code{data} based on the pre-processing
##' chosen for the PCA model
##' @param post unpre-process the final data (add the center back etc
##' to get the final estimate)
##' @param ... Not used
##' @return A matrix representing the fitted data
##' @keywords multivariate
##' @aliases fitted.pcaRes 
##' @examples
##' pc <- pca(iris[,1:4], nPcs=4, center=TRUE, scale="uv")
##' sum( (fitted(pc) - iris[,1:4])^2 )
##' @export
##' @author Henning Redestig
fitted.pcaRes <- function(object, data=NULL, nPcs=nP(object),
                          pre=TRUE, post=TRUE, ...) {
  if(method(object) %in% listPcaMethods("nonlinear")) {
    if(is.null(data) & is.null(completeObs(object)))
      stop("completeObs slot is empty -- provide the training data")
    if(is.null(data) & !is.null(completeObs(object)))
      data <- completeObs(object)
    if(is.null(data))
      stop("nlpca requires original data to be provide")
    if(pre)
      data <- prep(data, scl(object), center(object))
    recData <-
      errorHierarchic(object@network, t(scores(object)), t(data))$out[,,nPcs]
    recData <- t(recData / object@network@scalingFactor)
  }
  else  {
    if(!is.null(data)) {
      if(pre)
        data <- prep(data, scl(object), center(object))
      tt <- data %*% loadings(object)[,1:nPcs, drop=FALSE]
    }
    if(is.null(data))
      tt <- scores(object)[,1:nPcs, drop=FALSE]
    recData <- tcrossprod(tt, loadings(object)[,1:nPcs, drop=FALSE])
  }
  if(post)
    recData <- prep(recData, scl(object), center(object), reverse=TRUE)
  return(recData)
}
##' Fitted PCA data.
##' @seealso \code{\link{fitted.pcaRes}}
##' @importFrom stats,fitted
##' @exportMethod fitted
##' @aliases fitted,pcaRes-method
setMethod("fitted", "pcaRes", fitted.pcaRes)

#' Plot the computed diagnostics of PCA model to get an idea of their
##' importance. Note though that the standard screeplot shows the
##' standard deviations for the PCs this method shows the R2 values
##' which empirically shows the importance of the P's and is thus
##' applicable for any PCA method rather than just SVD based PCA.
##'
##' If cross-validation was done for the PCA the plot will also show
##' the CV based statistics.  A common rule-of-thumb for determining
##' the optimal number of PCs is the PC where the CV diagnostic is at
##' its maximum but not very far from \eqn{R^2}. 
##' @title Plot diagnostics (screeplot) 
##' @param x \code{pcaRes} The pcaRes object.
##' @param y not used
##' @param main title of the plot
##' @param col Colors of the bars
##' @param ... further arguments to barplot
##' @return None, used for side effect.
##' @seealso \link{screeplot}
##' @examples
##' data(metaboliteData)
##' pc <- pca(t(metaboliteData), nPcs=5, cv="q2", scale="uv")
##' plot(pc)
##' @export
##' @author Henning Redestig
plot.pcaRes <- function(x, y=NULL, main=deparse(substitute(object)),
                        col=gray(c(0.9, 0.5)), ...) {
  ## the deparse(subsitute(object)) later fails otherwise
  main <- main
  if(!is.null(cvstat(x))) {
    cvs <- cvstat(x)
    if(length(cvs) != nP(x))
      cvs <- c(cvs, rep(NA, nP(x) - length(cvs)))
    xx <- rbind(R2cum(x), cvs)
    barplot(xx, beside=TRUE, ylim=c(0,1.1), col=col, main=main,
            names.arg=paste("PC", 1:nP(x), sep=""), ...)
    legend(x="topleft", fill=col,
           legend=c(expression(R^2), expression(Q^2)))
  } else 
    barplot(R2cum(x), ylim=c(0,1.1), ylab=expression(R^2), main=main,
            names.arg=paste("PC", 1:nP(x), sep=""), col=col[1], ...)
}

##' Plot the R2 of the principal components to get an idea of their
##' importance. Note though that the standard screeplot shows the
##' standard deviations for the PC's this method shows the R2 values
##' which empirically shows the importance of the PC's and is thus
##' applicable for any PCA method rather than just SVD based PCA.
##' @note This method is deprecated in favor of plot.pcaRes which does
##' (almost) the same thing but with a better name.
##' @title R2 plot (screeplot) for PCA
##' @param object \code{pcaRes} The pcaRes object.
##' @param nPcs \code{numeric} The number of PC's to consider.
##' @param type \code{character} Barplot or line plot
##' @param main \code{character} The main label of the plot
##' @param ... Passed on to \code{screeplot}
##' @return None, used for side effect.
##' @seealso \link{screeplot}
##' @keywords multivariate
##' @export 
##' @author Henning Redestig
plotR2 <- function(object, nPcs=nP(object), type=c("barplot", "lines"),
                   main=deparse(substitute(object)), ...) {
  ## ditch from version 2.0.0
  .Deprecated("plot.pcaRes", "pcaMethods")
  main <- main
  newx <- list(sdev=object@R2)
  screeplot(newx, nPcs, type, main,...)
}

##' A common way of visualizing two principal components
##'
##' This method is meant to be used as a quick way to visualize
##' results, if you want a more specific plot you probably want to
##' get the scores, loadings with \code{scores(object)},
##' \code{loadings(object)} and then design your own plotting method.
##' @title Side by side scores and loadings plot
##' @usage slplot(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
##' sl="def", ll="def", hotelling=0.95, rug=TRUE, sub=NULL,...)
##' @param object a pcaRes object
##' @param pcs which two pcs to plot
##' @param scoresLoadings Which should be shown scores and or loadings
##' @param sl labels to plot in the scores plot
##' @param ll labels to plot in the loadings plot
##' @param hotelling confidence interval for ellipse in the score plot
##' @param rug logical, rug x axis in score plot or not
##' @param sub Subtitle, defaults to annotate with amount of explained
##' variance.
##' @param ... Further arguments to plot functions. Prefix arguments
##' to \code{par()} with 's' for the scores plot and 'l' for the
##' loadings plot. I.e. cex become scex for setting character
##' expansion in the score plot and lcex for the loadings plot.
##' @return None, used for side effect.
##' @note Uses layout instead of par to provide side-by-side so it
##' works with Sweave (but can not be combined with
##' \code{par(mfrow=..))}
##' @author Henning Redestig
##' @seealso \code{\link{pca}}, \code{\link{biplot}}
##' @exportMethod slplot
##' @aliases slplot slplot,pcaRes-method
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4], scale="uv")
##' slplot(pcIr, sl=NULL, spch=5)
##' slplot(pcIr, sl=NULL, lcex=1.3, scol=as.integer(iris[,5]))
##' @keywords multivariate
setMethod("slplot", "pcaRes",
function(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
         sl=rownames(scores(object)),
         ll=rownames(loadings(object)), hotelling=0.95, rug=FALSE,
         sub=NULL,...) {

  opar <- par(no.readonly=TRUE)

  cl <- match.call()
  mainArgs <- c(1,match(c("ll", "sl", "scoresLoadings", "sub"),
                        names(cl), 0))
  scoreArgs <- grep("^s", names(cl)[-mainArgs])
  loadingArgs <- grep("^l", names(cl)[-mainArgs])

  if(!is.null(ll) & length(ll) != nVar(object))
    stop("Loading labels do not match the object dimensions")
  if(!is.null(sl) & length(sl) != nObs(object))
    stop("Score labels do not match the object dimensions")
  if(is.null(sl))
    sl <- NA
  if(is.null(ll))
    ll <- NA
  
  ## no loadings for non-linear pca
  if(method(object) %in% listPcaMethods("nonlinear") && scoresLoadings[2])
    scoresLoadings[2] <- FALSE
  
  if(length(pcs) > 2)
    plotPcs(object, pcs, scoresLoadings=scoresLoadings,...)
  else {
    if(is.null(sub))
      sub <- paste(sprintf("%.2f", R2cum(object)[max(pcs)]
                           * 100),
                   "% of the variance explained", sep="")

    if(sum(scoresLoadings) == 2)
      layout(matrix(c(1,2), 1, 2, TRUE), respect=matrix(c(1,1), 1, 2))
    ## exception plot if one dimensional
    if (length(pcs) == 1 | nP(object) == 1) {
      pcs <- 1
      
      ## score plot
      if(scoresLoadings[1]) {
        newCall <- call("barplot",
                        height=scores(object)[,pcs],
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
                        height=loadings(object)[,pcs],
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
                       x=scores(object)[,pcs],
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
                         x=scores(object)[,pcs], labels=sl)
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
        rug(scores(object)[,1])
      abline(h=0, v=0)
      if(!is.null(hotelling)) {
        A <- length(pcs)
        el <- simpleEllipse(scores(object)[,pcs[1]],
                            scores(object)[,pcs[2]], alfa=hotelling)
        lines(el)
      }
    }

    ## the loading plot
    if(scoresLoadings[2]) {
      ## setup plot
      plotCall <- call("plot",
                       x=loadings(object)[,pcs],
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
                         x=loadings(object)[,pcs], labels=ll)
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

