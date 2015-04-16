##' The leverages of PCA model indicate how much influence each
##' observation has on the PCA model. Observations with high leverage
##' has caused the principal components to rotate towards them. It can
##' be used to extract both "unimportant" observations as well as
##' picking potential outliers.
##'
##' Defined as \eqn{Tr(T(T'T)^{-1}T')}{Tr(T(T'T)^(-1)T')}
##' @title Extract leverages of a PCA model
##' @param object a \code{pcaRes} object
##' @return The observation leverages as a numeric vector
##' @references Introduction to Multi- and Megavariate Data Analysis
##' using Projection Methods (PCA and PLS), L. Eriksson, E. Johansson,
##' N. Kettaneh-Wold and S. Wold, Umetrics 1999, p. 466
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4])
##' ## versicolor has the lowest leverage
##' with(iris, plot(leverage(pcIr)~Species))
##' @keywords multivariate
##' @aliases leverage leverage,pcaRes-method
##' @author Henning Redestig
setMethod("leverage", "pcaRes",
          function(object) {
            diag(scores(object) %*%
                 solve(crossprod(scores(object))) %*% t(scores(object)))
          })

##' Distance to the model of X-space.
##'
##' Measures how well described the observations are, i.e. how well
##' they fit in the mode. High DModX indicate a poor fit. Defined as:
##'
##' \eqn{\frac{\sqrt{\frac{SSE_i}{K-A}}}{\sqrt{\frac{SSE}{(N-A-A_0)(K-A)}}}}
##'
##' For observation \eqn{i}, in a model with \eqn{A} components,
##' \eqn{K} variables and \eqn{N} obserations. SSE is the squared sum
##' of the residuals. \eqn{A_0} is 1 if model was centered and 0
##' otherwise. DModX is claimed to be approximately F-distributed and
##' can therefore be used to check if an observation is significantly
##' far away from the PCA model assuming normally distributed data.
##'
##' Pass original data as an argument if the model was calculated with
##' \code{completeObs=FALSE}.
##' @title DModX
##' @usage DModX(object, dat, newdata=FALSE, type=c("normalized","absolute"), ...)
##' @param object a pcaRes object
##' @param dat the original data, taken from \code{completeObs} if
##' left missing.
##' @param newdata logical indicating if this data was part of the
##' training data or not. If it was, it is adjusted by a near one factor
##' \eqn{v=(N/ (N-A-A0))^-1}
##' @param type if absolute or normalized values should be
##' given. Normalized values are adjusted to the the total RSD of the
##' model.
##' @param ... Not used 
##' @return A vector with distances from observations to the PCA model
##' @aliases DModX DModX,pcaRes-method
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4])
##' with(iris, plot(DModX(pcIr)~Species))
##' @references Introduction to Multi- and Megavariate Data Analysis
##' using Projection Methods (PCA and PLS), L. Eriksson, E. Johansson,
##' N. Kettaneh-Wold and S. Wold, Umetrics 1999, p. 468
##' @author Henning Redestig
setMethod("DModX", "pcaRes",
          function(object, dat, newdata=FALSE, type=c("normalized","absolute"), ...) {
            type <- match.arg(type)
            if(missing(dat)) {
              if(!is.null(completeObs(object)))
                dat <- completeObs(object)
              else stop("missing data when calculating DModX")
            }
            A0 <- as.integer(centered(object))
            ny <- ifelse(newdata, 1,
                         sqrt(nObs(object) /
                              (nObs(object) - nP(object) - A0)))
            E2 <- resid(object, dat)^2
            s <- sqrt(rowSums(E2) / (nVar(object) - nP(object))) * ny
            if(type == "absolute")
              return(s)
            s0 <- sqrt(sum(E2) / ((nObs(object) - nP(object) - A0) *
                                  (nVar(object) - nP(object))))
            s / s0
          })

##' Get number of PCs
##' @param object pcaRes object
##' @param ...  not used
##' @return Number of PCs
##' @aliases nP nP,pcaRes-method
##' @usage nP(object, ...)
##' @author Henning Redestig
setMethod("nP", "pcaRes", function(object, ...) {
  if(is.null(object@nPcs) & !is.null(scores(object)))
    return(ncol(scores(object)))
  object@nPcs
})


##' Get cross-validation statistics (e.g. \eqn{Q^2}).
##' @param object pcaRes object
##' @param ...  not used
##' @return vector CV statistics 
##' @aliases cvstat cvstat,pcaRes-method
##' @usage cvstat(object, ...)
##' @author Henning Redestig
setMethod("cvstat", "pcaRes", function(object, ...) {
  object@cvstat
})

##' Get number of PCs.
##' @param object pcaRes object
##' @param ...  not used
##' @note Try to use \code{link{nP}} instead since \code{nPcs} tend to
##' clash with argument names.
##' @return Number of PCs
##' @usage nPcs(object, ...)
##' @aliases nPcs nPcs,pcaRes-method
##' @author Henning Redestig
setMethod("nPcs", "pcaRes", function(object, ...) {
  nP(object)
})

##' Get the number of observations used to build the PCA model.
##' @param object pcaRes object
##' @param ... Not used
##' @usage nObs(object, ...)
##' @aliases nObs nObs,pcaRes-method
##' @return Number of observations
##' @author Henning Redestig
setMethod("nObs", "pcaRes", function(object, ...) {
  object@nObs
})

##' Get the number of variables used to build the PCA model.
##' @param object pcaRes object
##' @param ... Not used
##' @usage nVar(object, ...)
##' @aliases nVar nVar,pcaRes-method
##' @return Number of variables
##' @author Henning Redestig
setMethod("nVar", "pcaRes", function(object, ...) {
  object@nVar
})

##' Check centering was part of the model
##' @param object pcaRes object
##' @param ... Not used
##' @usage centered(object, ...)
##' @aliases centered centered,pcaRes-method
##' @return TRUE if model was centered
##' @author Henning Redestig
setMethod("centered", "pcaRes", function(object, ...) {
  if(is.null(object@centered))
    return(FALSE)
  object@centered  
})

##' Get the centers of the original variables
##' @param object pcaRes object
##' @param ... Not used
##' @usage center(object, ...)
##' @aliases center center,pcaRes-method
##' @return Vector with the centers
##' @author Henning Redestig
setMethod("center", "pcaRes", function(object, ...) {
  object@center
})

setMethod("completeObs", "pcaRes", function(object, ...) {
  object@completeObs
})
##' Get the original data with missing values replaced with predicted
##' values.
##' @param object object to fetch complete data from
##' @param ... Not used
##' @usage completeObs(object, ...)
##' @aliases completeObs completeObs,nniRes-method
##' completeObs,pcaRes-method
##' @return Completed data (matrix)
##' @author Henning Redestig
setMethod("completeObs", "nniRes", function(object, ...) {
  object@completeObs
})

##' Get the used PCA method
##' @param object pcaRes object
##' @param ... Not used
##' @usage method(object, ...)
##' @aliases method method,pcaRes-method
##' @return The used pca method
##' @author Henning Redestig
setMethod("method", "pcaRes", function(object, ...) {
  object@method
})

setMethod("nmissing", "nniRes", function(object, ...) {
  sum(object@missing)
})
##' Missing values
##' @param object pcaRes object
##' @param ... Not used
##' @usage nmissing(object, ...)
##' @aliases nmissing nmissing,pcaRes-method nmissing,nniRes-method
##' @return Get the number of missing values
##' @author Henning Redestig
setMethod("nmissing", "pcaRes", function(object, ...) {
  sum(object@missing)
})

##' Get a matrix with indicating the elements that were missing in the
##' input data. Convenient for estimating imputation performance.
##' @param object pcaRes object
##' @param ... Not used
##' @usage wasna(object, ...)
##' @aliases wasna wasna,pcaRes-method
##' @return A matrix with logicals
##' @examples
##' data(metaboliteData)
##' data(metaboliteDataComplete)
##' result <- pca(metaboliteData, nPcs=2)
##' plot(completeObs(result)[wasna(result)], metaboliteDataComplete[wasna(result)])
##' @author Henning Redestig
setMethod("wasna", "pcaRes", function(object, ...) {
  object@missing
})

##' Get the standard deviations of the scores (indicates their
##' relevance)
##' @param object pcaRes object
##' @param ... Not used
##' @usage sDev(object, ...)
##' @aliases sDev sDev,pcaRes-method
##' @return Standard devations of the scores
##' @author Henning Redestig
setMethod("sDev", "pcaRes", function(object, ...) {
  object@sDev
})

##' Check if scaling was part of the PCA model
##' @param object pcaRes object
##' @param ... Not used
##' @usage scaled(object, ...)
##' @aliases scaled scaled,pcaRes-method
##' @return TRUE if scaling was part of the PCA model
##' @author Henning Redestig
setMethod("scaled", "pcaRes", function(object, ...) {
  if(is.null(object@scaled))
    return(FALSE)
  object@scaled != "none"
})

##' Get the scales (e.g. standard deviations) of the original
##' variables
##' @param object pcaRes object
##' @param ... Not used
##' @usage scl(object, ...)
##' @aliases scl scl,pcaRes-method
##' @return Vector with the scales
##' @seealso \code{\link{prep}}
##' @author Henning Redestig
setMethod("scl", "pcaRes", function(object, ...) {
  object@scale 
})

##' Cumulative R2 is the total ratio of variance that is being
##' explained by the model
##' @param object a \code{pcaRes} model
##' @param ...  Not used
##' @return Get the cumulative R2
##' @aliases R2cum R2cum,pcaRes-method
##' @author Henning Redestig
setMethod("R2cum", "pcaRes", function(object, ...) {
  object@R2cum
})

##' Get scores from a pcaRes object
##' @param object a pcaRes object
##' @param ... not used
##' @return The scores as a matrix
##' @export
##' @author Henning Redestig
##' @method scores pcaRes
scores.pcaRes <- function(object, ...) object@scores
##' Get scores from a pcaRes object
##' @param object a pcaRes object
##' @param ... not used
##' @return The scores as a matrix
##' @seealso \code{\link{scores.pcaRes}}
##' @aliases scores scores,pcaRes-method
##' @author Henning Redestig
setMethod("scores", "pcaRes", scores.pcaRes)

##' Get loadings from a pcaRes object
##' @param object a pcaRes object
##' @param ... not used
##' @return The loadings as a matrix
##' @export
##' @author Henning Redestig
##' @method loadings pcaRes
loadings.pcaRes <- function(object, ...) object@loadings
##' Get loadings from a pcaRes object
##' @param object a pcaRes object
##' @param ...  not used
##' @return The loadings as a matrix
##' @seealso \code{\link{loadings.pcaRes}}
##' @author Henning Redestig
##' @aliases loadings,pcaRes-method
setMethod("loadings", "pcaRes", loadings.pcaRes)
##' Crude way to unmask the function with the same name from
##' \code{stats}
##' @param object any object
##' @param ...  not used
##' @return The loadings
##' @author Henning Redestig
##' @aliases loadings loadings,ANY-method
setMethod("loadings", "ANY", function(object,...) {
  stats::loadings(object)
})


##' Dimensions of a PCA model
##' @param x a pcaRes object
##' @return Get the dimensions of this PCA model
##' @method dim pcaRes
##' @export 
##' @author Henning Redestig
dim.pcaRes <- function(x)  {
  res <-  c(nObs(x), nVar(x), nP(x))
  names(res) <- c("nObs", "nVar", "nPcs")
  res
}

##' Print basic information about pcaRes object
##' @title Print/Show for pcaRes
##' @param x a pcaRes object
##' @param ...  not used
##' @return nothing, used for its side effect
##' @name show-methods
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
##' @aliases print,pcaRes-method print,nniRes-method
##' @name show-methods
setMethod("print", "pcaRes", showPcaRes)
## @importFrom methods show
##' @aliases show,pcaRes-method show,nniRes-method
##' @param object the object to print information about
##' @name show-methods
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
##' @method biplot pcaRes
##' @export
##' @examples
##' data(iris)
##' pcIr <- pca(iris[,1:4])
##' biplot(pcIr)
##' @seealso \code{prcomp}, \code{pca}, \code{princomp}
##' @author Kevin Wright, Adapted from \code{biplot.prcomp}
##' @keywords multivariate
##' @name biplot-methods
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
##' @aliases biplot,pcaRes-method
## @importFrom stats biplot
##' @name biplot-methods
setMethod("biplot", "pcaRes", biplot.pcaRes)

##' Flexible calculation of R2 goodness of fit.
##' @title R2 goodness of fit
##' @param object a PCA model object
##' @param direction choose between calculating R2 per variable, per
##' observation or for the entire data with 'variables',
##' 'observations' or 'complete'.
##' @param data the data used to fit the model
##' @param pcs the number of PCs to use to calculate R2
##' @aliases R2VX R2VX,pcaRes-method
##' @examples
##' R2VX(pca(iris))
##' @return A vector with R2 values
##' @author Henning Redestig
setMethod('R2VX', 'pcaRes',
          function(object,
                   direction=c('variables', 'observations', 'complete'),
                   data=completeObs(object), pcs=nP(object)) {
  direction <- match.arg(direction)
  if(is.null(data)) 
    stop('missing input when calculating R2')
  if(any(is.na(data))) 
    stop('missing values not allowed for calculating R2')
  dat <- prep(data, scale=scl(object), center=center(object))
  xhat <- resid(object, pcs=pcs)
  switch(direction, variables={
    1 - colSums(xhat^2) / colSums(dat^2)
  }, observations={
    1 - rowSums(xhat^2) / rowSums(dat^2)
  }, complete={
    1 - sum(xhat^2) / sum(dat^2)
  })
})

setAs('pcaRes', 'data.frame', function(from) {
  tt <- scores(from)
  pp <- loadings(from)
  if(is.null(rownames(tt)))
    rownames(tt) <- 1:nrow(tt)
  if(is.null(rownames(pp)))
    rownames(pp) <- 1:nrow(pp)
  dfs <- as.data.frame(tt)
  dfs$names <- rownames(tt)
  dfs$type <- 'scores'
  dfl <- as.data.frame(pp)
  dfl$names <- rownames(pp)
  dfl$type <- 'loadings'
  rownames(dfl) <- rownames(dfs) <- NULL
  rbind(dfl, dfs)
})

##' Print a brief description of the PCA model
##' @title Summary of PCA model
##' @param object a pcaRes object
##' @param ... Not used
##' @return Nothing, used for side-effect
##' @aliases summary summary.pcaRes summary,pcaRes-method
##' @author Henning Redestig
##' @export
##' @name summary
##' @method summary pcaRes
summary.pcaRes <- function(object, ...){
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
}
setMethod("summary", "pcaRes", summary.pcaRes)

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
##' @method predict pcaRes 
##' @keywords multivariate
##' @examples
##' data(iris)
##' hidden <- sample(nrow(iris), 50)
##' pcIr <- pca(iris[-hidden,1:4])
##' pcFull <- pca(iris[,1:4])
##' irisHat <- predict(pcIr, iris[hidden,1:4])
##' cor(irisHat$scores[,1], scores(pcFull)[hidden,1])
##' @export
##' @name predict-methods
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
## @importFrom stats predict
##' @name predict-methods
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
##' @method residuals pcaRes
##' @keywords multivariate
##' @export
##' @name rediduals-methods
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
##' @aliases residuals,pcaRes-method
##' @name rediduals-methods
setMethod("residuals", "pcaRes", residuals.pcaRes)
##' @name rediduals-methods
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
##' @method fitted pcaRes
##' @examples
##' pc <- pca(iris[,1:4], nPcs=4, center=TRUE, scale="uv")
##' sum( (fitted(pc) - iris[,1:4])^2 )
##' @export
##' @name fitted-methods
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
## @importFrom stats fitted
##' @name fitted-methods
##' @aliases fitted,pcaRes-method
setMethod("fitted", "pcaRes", fitted.pcaRes)

##' Plot the computed diagnostics of PCA model to get an idea of their
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
##' @method plot pcaRes
##' @aliases plot.pcaRes plot,pcaRes-method
##' @export
##' @author Henning Redestig
plot.pcaRes <- function(x, y=NULL, main=deparse(substitute(object)),
                        col=gray(c(0.9, 0.5)), ...) {
  y <- NULL
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
setMethod("plot", "pcaRes", plot.pcaRes)

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
