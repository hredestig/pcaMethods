setGeneric("slplot", function(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
                              sl="def", ll="def",
                              hotelling=0.95, rug=TRUE, sub=NULL,...) standardGeneric("slplot"))

setGeneric("vector2matrices", function(object, ...) standardGeneric("vector2matrices"))

##setGeneric("biplot", function(x, ...) standardGeneric("biplot"))
##setGeneric("biplot", function(x, choices=1:2, scale=1, pc.biplot=FALSE, ...) standardGeneric("biplot"))
