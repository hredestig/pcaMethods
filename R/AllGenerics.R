setGeneric("vector2matrices", function(object, ...)
           standardGeneric("vector2matrices"))
##' @exportMethod leverage
setGeneric("leverage", function(object, ...) standardGeneric("leverage"))
##' @exportMethod DModX
setGeneric("DModX", function(object, dat, newdata=FALSE, type=c("normalized","absolute"), ...)
    standardGeneric("DModX"))
##' @exportMethod nP
setGeneric("nP", function(object, ...) standardGeneric("nP"))
##' @exportMethod cvstat
setGeneric("cvstat", function(object, ...) standardGeneric("cvstat"))
##' @exportMethod nPcs
setGeneric("nPcs", function(object, ...) standardGeneric("nPcs"))
##' @exportMethod nObs 
setGeneric("nObs", function(object, ...) standardGeneric("nObs"))
##' @exportMethod nVar 
setGeneric("nVar", function(object, ...) standardGeneric("nVar"))
##' @exportMethod centered 
setGeneric("centered", function(object, ...) standardGeneric("centered"))
##' @exportMethod center
setGeneric("center", function(object, ...) standardGeneric("center"))
##' @exportMethod completeObs
setGeneric("completeObs", function(object, ...) standardGeneric("completeObs"))
##' @exportMethod method
setGeneric("method", function(object, ...) standardGeneric("method"))
##' @exportMethod nmissing
setGeneric("nmissing", function(object, ...) standardGeneric("nmissing"))
##' @exportMethod wasna
setGeneric("wasna", function(object, ...) standardGeneric("wasna"))
##' @exportMethod sDev
setGeneric("sDev", function(object, ...) standardGeneric("sDev"))
##' @exportMethod scaled
setGeneric("scaled", function(object, ...) standardGeneric("scaled"))
##' @exportMethod scl
setGeneric("scl", function(object, ...) standardGeneric("scl"))
##' @exportMethod R2cum 
setGeneric("R2cum", function(object, ...) standardGeneric("R2cum"))
##' @exportMethod slplot
setGeneric("slplot", function(object, pcs=c(1,2),
                              scoresLoadings=c(TRUE, TRUE),
                              sl="def", ll="def",
                              hotelling=0.95, rug=TRUE, sub=NULL,...)
           standardGeneric("slplot"))
##' @exportMethod scores
setGeneric("scores", function(object, ...) standardGeneric("scores"))
##' @exportMethod loadings
setGeneric("loadings", function(object, ...) standardGeneric("loadings"))
##' @exportMethod R2VX
setGeneric("R2VX", function(object, ...) standardGeneric("R2VX"))
## @exportMethod prep
#setGeneric("prep", function(object, ...) standardGeneric("prep"))
