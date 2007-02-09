## class representation of the NLPCA neural net
setClass("nlpcaNet",
         representation(net="matrix",
                        hierarchic="list",
                        fct="character",
                        fkt="character",
                        weightDecay="numeric",
                        featureSorting="logical",
                        dataDist="matrix",
                        inverse="logical",
                        fCount="integer",
                        componentLayer="integer",
                        error="function",
                        gradient="function",
                        weights="list",
                        maxIter="integer",
                        scalingFactor="numeric"),
         prototype(net=rbind(c(4,6,2,6,4)),
                   hierarchic=list(var=rbind(c(1,1,0.01)), layer=3, idx=rbind(c(1,1,0),c(0,1,1))),
                   fct=c("linr", "tanh", "linr", "tanh", "linr"),
                   fkt=c("tanh", "linr", "tanh", "linr"),
                   weightDecay=0.001,
                   featureSorting=TRUE,
                   inverse=FALSE,
                   dataDist=NULL,
                   fCount=as.integer(0),
                   componentLayer=as.integer(3),
                   error=errorHierarchic,
                   gradient=derrorHierarchic,
                   weights=NULL,
                   maxIter=as.integer(1200),
                   scalingFactor=NULL))

##
## pcaRes is used by all PCA-based methods
##
setClass("pcaRes",
         representation(completeObs="matrix",
			scores="matrix",
                        loadings="matrix",
                        R2cum="numeric",
                        R2="numeric",
                        sDev="numeric",
                        nObs="numeric",
                        nVar="numeric",
                        centered="logical",
                        center="numeric",
                        varLimit="numeric",
                        nPcs="numeric",
                        method="character",
                        missing="numeric",
                        network="nlpcaNet"),
         prototype(completeObs=NULL,
		   scores=NULL,
                   loadings=NULL,
                   R2cum=NULL,
                   R2=NULL,
                   sDev=NULL,
                   nObs=NULL,
                   nVar=NULL,
                   centered=NULL,
                   center=NULL,
                   varLimit=NULL,
                   nPcs=NULL,
                   method=NULL,
                   missing=NULL,
                   network=NULL))
                   
setAs("NULL", "pcaRes",
      function(from, to){
        new(to)
      })


##
## clusterRes is used by all cluster based methods
##
setClass("nniRes",
         representation(completeObs="matrix",
                        nObs="numeric",
                        nVar="numeric",
                        centered="logical",
                        center="numeric",
                        k="numeric",
                        method="character",
                        correlation="character",
                        missing="numeric"),
         prototype(completeObs=NULL,
                   nObs=NULL,
                   nVar=NULL,
                   centered=NULL,
                   center=NULL,
                   k=NULL,
                   method=NULL,
                   correlation=NULL,
                   missing=NULL))

setAs("NULL", "nniRes",
      function(from, to) {
          new(to)
      })



weightsAccount <- function(w) {
  list(
       set = function(newWeights) {
         if(!inherits(newWeights, "matrix"))
           stop("The weights must inherit from matrix")
          w <<- newWeights
       },
       current = function() {
         w
       }
       )
}

setAs("NULL", "nlpcaNet",
      function(from, to){
        new(to)
      })

