setGeneric("slplot", function(object, pcs=c(1,2),
                              scoresLoadings=c(TRUE, TRUE),
                              sl="def", ll="def",
                              hotelling=0.95, rug=TRUE, sub=NULL,...)
           standardGeneric("slplot"))
setGeneric("vector2matrices", function(object, ...)
           standardGeneric("vector2matrices"))

setGeneric("leverage", function(object, ...) standardGeneric("leverage"))
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
##' @exportMethod leverage
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

setGeneric("DModX", function(object, dat, ...) standardGeneric("DModX"))
##' Distance to the model of X-space.
##'
##' Measures how well described the observations are, i.e. how well
##' they fit in the mode. High DModX indicate a poor fit. Defined as:
##'
##' \eqn{\frac{\sqrt{\frac{SSE_i}{K-A}}}{\sqrt{\frac{SSE}{(N-A-A_0)(K-A)}}}}
##'
##' For observation \eqn{i}, in a model with \eqn{A} components,
##' \eqn{K} variables and \eqn{N} obserations. SSE is the squared sum
##' of the residuals. \eqn{A_0} is 0 if model was centered and 1
##' otherwise. DModX is claimed to be approximately F-distributed and
##' can therefore be used to check if an observation is significantly
##' far away from the PCA model assuming normally distributed data.
##'
##' Pass original data as an argument if the model was calculated with
##' \code{completeObs=FALSE}.
##' @title DModX
##' @usage DModX(object, dat, ...)
##' @param object a pcaRes object
##' @param dat the original data, taken from \code{completeObs} if
##' left missing.
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
##' @exportMethod DModX
##' @author Henning Redestig
setMethod("DModX", "pcaRes",
          function(object, dat, ...) {
            newdata <- TRUE
            if(missing(dat)) {
              newdata <- FALSE
              if(!is.null(completeObs(object))) dat <- completeObs(object)
              else stop("missing data when calculating DModX")
            }
            newE2 <- resid(object, dat)^2
            modelE2 <- resid(object, completeObs(object))^2
            sEk2 <- sqrt(newE2  / (nVar(object) - nP(object)))
            ny <- sqrt(nObs(object) /
                       (nObs(object) - nP(object) - as.integer(centered(object))))
            s0 <- sqrt(sum(modelE2) /
                       ((nObs(object) - nP(object) - as.integer(centered(object))) *
                        (nVar(object) - nP(object))))
            if(newdata)
              ny <- 1
            (rowSums(sEk2, na.rm=TRUE) * ny) / s0
          })

setGeneric("nP", function(object, ...) standardGeneric("nP"))
##' Get number of PCs
##' @param object pcaRes object
##' @param ...  not used
##' @return Number of PCs
##' @exportMethod nP
##' @aliases nP nP,pcaRes-method
##' @usage nP(object, ...)
##' @author Henning Redestig
setMethod("nP", "pcaRes", function(object, ...) {
  if(is.null(object@nPcs) & !is.null(scores(object)))
    return(ncol(scores(object)))
  object@nPcs
})


setGeneric("cvstat", function(object, ...) standardGeneric("cvstat"))
##' Get cross-validation statistics (e.g. \eqn{Q^2}).
##' @param object pcaRes object
##' @param ...  not used
##' @return vector CV statistics 
##' @exportMethod cvstat
##' @aliases cvstat cvstat,pcaRes-method
##' @usage cvstat(object, ...)
##' @author Henning Redestig
setMethod("cvstat", "pcaRes", function(object, ...) {
  object@cvstat
})

setGeneric("nPcs", function(object, ...) standardGeneric("nPcs"))
##' Get number of PCs.
##' @param object pcaRes object
##' @param ...  not used
##' @note Try to use \code{link{nP}} instead since \code{nPcs} tend to
##' clash with argument names.
##' @return Number of PCs
##' @exportMethod nPcs
##' @usage nPcs(object, ...)
##' @aliases nPcs nPcs,pcaRes-method
##' @author Henning Redestig
setMethod("nPcs", "pcaRes", function(object, ...) {
  nP(object)
})

setGeneric("nObs", function(object, ...) standardGeneric("nObs"))
##' Get the number of observations used to build the PCA model.
##' @param object pcaRes object
##' @param ...
##' @usage nObs(object, ...)
##' @aliases nObs nObs,pcaRes-method
##' @return Number of observations
##' @exportMethod nObs 
##' @author Henning Redestig
setMethod("nObs", "pcaRes", function(object, ...) {
  object@nObs
})

setGeneric("nVar", function(object, ...) standardGeneric("nVar"))
##' Get the number of variables used to build the PCA model.
##' @param object pcaRes object
##' @param ...
##' @usage nVar(object, ...)
##' @aliases nVar nVar,pcaRes-method
##' @return Number of variables
##' @exportMethod nVar 
##' @author Henning Redestig
setMethod("nVar", "pcaRes", function(object, ...) {
  object@nVar
})

setGeneric("centered", function(object, ...) standardGeneric("centered"))
##' Check centering was part of the model
##' @param object pcaRes object
##' @param ... Not used
##' @usage centered(object, ...)
##' @aliases centered centered,pcaRes-method
##' @return TRUE if model was centered
##' @exportMethod centered 
##' @author Henning Redestig
setMethod("centered", "pcaRes", function(object, ...) {
  if(is.null(object@centered))
    return(FALSE)
  object@centered  
})

setGeneric("center", function(object, ...) standardGeneric("center"))
##' Get the centers of the original variables
##' @param object pcaRes object
##' @param ... Not used
##' @usage center(object, ...)
##' @aliases center center,pcaRes-method
##' @return Vector with the centers
##' @exportMethod center
##' @author Henning Redestig
setMethod("center", "pcaRes", function(object, ...) {
  object@center
})

setGeneric("completeObs", function(object, ...) standardGeneric("completeObs"))
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
##' @exportMethod completeObs
##' @author Henning Redestig
setMethod("completeObs", "nniRes", function(object, ...) {
  object@completeObs
})

setGeneric("method", function(object, ...) standardGeneric("method"))
##' Get the used PCA method
##' @param object pcaRes object
##' @param ... Not used
##' @usage method(object, ...)
##' @aliases method method,pcaRes-method
##' @return The used pca method
##' @exportMethod method
##' @author Henning Redestig
setMethod("method", "pcaRes", function(object, ...) {
  object@method
})

setGeneric("nmissing", function(object, ...) standardGeneric("nmissing"))
setMethod("nmissing", "nniRes", function(object, ...) {
  sum(object@missing)
})
##' Missing values
##' @param object pcaRes object
##' @param ... Not used
##' @usage nmissing(object, ...)
##' @aliases nmissing nmissing,pcaRes-method nmissing,nniRes-method
##' @return Get the number of missing values
##' @exportMethod nmissing
##' @author Henning Redestig
setMethod("nmissing", "pcaRes", function(object, ...) {
  sum(object@missing)
})

setGeneric("wasna", function(object, ...) standardGeneric("wasna"))
##' Get a matrix with indicating the elements that were missing in the
##' input data. Convenient for estimating imputation performance.
##' @param object pcaRes object
##' @param ... Not used
##' @usage wasna(object, ...)
##' @aliases wasna wasna,pcaRes-method
##' @return A matrix with logicals
##' @exportMethod wasna
##' @examples
##' data(metaboliteData)
##' data(metaboliteDataComplete)
##' result <- pca(metaboliteData, nPcs=2)
##' plot(completeObs(result)[wasna(result)], metaboliteDataComplete[wasna(result)])
##' @author Henning Redestig
setMethod("wasna", "pcaRes", function(object, ...) {
  object@missing
})

setGeneric("sDev", function(object, ...) standardGeneric("sDev"))
##' Get the standard deviations of the scores (indicates their
##' relevance)
##' @param object pcaRes object
##' @param ... Not used
##' @usage sDev(object, ...)
##' @aliases sDev sDev,pcaRes-method
##' @return Standard devations of the scores
##' @exportMethod sDev
##' @author Henning Redestig
setMethod("sDev", "pcaRes", function(object, ...) {
  object@sDev
})

setGeneric("scaled", function(object, ...) standardGeneric("scaled"))
##' Check if scaling was part of the PCA model
##' @param object pcaRes object
##' @param ... Not used
##' @usage scaled(object, ...)
##' @aliases scaled scaled,pcaRes-method
##' @return TRUE if scaling was part of the PCA model
##' @exportMethod scaled
##' @author Henning Redestig
setMethod("scaled", "pcaRes", function(object, ...) {
  if(is.null(object@scaled))
    return(FALSE)
  object@scaled != "none"
})

setGeneric("scl", function(object, ...) standardGeneric("scl"))
##' Get the scales (e.g. standard deviations) of the original
##' variables
##' @param object pcaRes object
##' @param ... Not used
##' @usage scl(object, ...)
##' @aliases scl scl,pcaRes-method
##' @return Vector with the scales
##' @seealso \code{\link{prep}}
##' @exportMethod scl
##' @author Henning Redestig
setMethod("scl", "pcaRes", function(object, ...) {
  object@scale 
})

setGeneric("R2cum", function(object, ...) standardGeneric("R2cum"))
##' Cumulative R2 is the total ratio of variance that is being
##' explained by the model
##' @param object a \code{pcaRes} model
##' @param ...  Not used
##' @return Get the cumulative R2
##' @aliases R2cum R2cum,pcaRes-method
##' @exportMethod R2cum 
##' @author Henning Redestig
setMethod("R2cum", "pcaRes", function(object, ...) {
  object@R2cum
})

##' Get the scores from a PCA model
##' @param object a pcaRes object
##' @param ... not used
##' @return The scores as a matrix
##' @export 
##' @author Henning Redestig
scores.pcaRes <- function(object,...) 
  object@scores

##' Get the loadings from a PCA model
##' @param object a pcaRes object
##' @param ...  not used
##' @return The loadings as a matrix
##' @export 
##' @author Henning Redestig
loadings.pcaRes <- function(object,...) 
  object@loadings

##' Dimensions of a PCA model
##' @param x a pcaRes object
##' @return Get the dimensions of this PCA model
##' @export 
##' @author Henning Redestig
dim.pcaRes <- function(x)  {
  res <-  c(nObs(x), nVar(x), nP(x))
  names(res) <- c("nObs", "nVar", "nPcs")
  res
}
