setGeneric("slplot", function(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
                              sl="def", ll="def",
                              hotelling=0.95, rug=TRUE, sub=NULL,...)
           standardGeneric("slplot"))

setGeneric("vector2matrices", function(object, ...) standardGeneric("vector2matrices"))
setGeneric("leverage", function(object, ...) standardGeneric("leverage"))

## getters
setGeneric("nPcs", function(object, ...) standardGeneric("nPcs"))
setGeneric("nObs", function(object, ...) standardGeneric("nObs"))
setGeneric("nVar", function(object, ...) standardGeneric("nVar"))
setGeneric("completeObs", function(object, ...) standardGeneric("completeObs"))
setGeneric("centered", function(object, ...) standardGeneric("centered"))
setGeneric("method", function(object, ...) standardGeneric("method"))
setGeneric("sDev", function(object, ...) standardGeneric("sDev"))
setGeneric("scaled", function(object, ...) standardGeneric("scaled"))
setGeneric("center", function(object, ...) standardGeneric("center"))


##setGeneric("biplot", function(x, ...) standardGeneric("biplot"))
##setGeneric("biplot", function(x, choices=1:2, scale=1, pc.biplot=FALSE, ...) standardGeneric("biplot"))
