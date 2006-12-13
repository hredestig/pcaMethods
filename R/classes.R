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
                        subset="numeric",
                        varLimit="numeric",
                        scaled="character",
                        nPcs="numeric",
                        method="character",
                        missing="numeric"),
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
                   subset=NULL,
                   varLimit=NULL,
                   scaled=NULL,
                   nPcs=NULL,
                   method=NULL,
                   missing=NULL))
                   
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
                        subset="numeric",
                        scaled="character",
                        k="numeric",
                        method="character",
                        correlation="character",
                        missing="numeric"),
         prototype(completeObs=NULL,
                   nObs=NULL,
                   nVar=NULL,
                   centered=NULL,
                   center=NULL,
                   subset=NULL,
                   scaled=NULL,
                   k=NULL,
                   method=NULL,
                   correlation=NULL,
                   missing=NULL))

setAs("NULL", "nniRes",
      function(from, to) {
          new(to)
      })

