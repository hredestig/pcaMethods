##' @include errorHierarchic.R
##' @include derrorHierarchic.R
NULL

##' This is a class representation of a non-linear PCA neural
##' network. The \code{nlpcaNet} class is not meant for user-level
##' usage.
##'
##' Creating Objects
##' 
##' \code{new("nlpcaNet", net=[the network structure],
##' hierarchic=[hierarchic design],
##' fct=[the functions at each layer], fkt=[the functions used for
##' forward propagation], weightDecay=[incremental decrease of weight
##' changes over iterations (between 0 and 1)], featureSorting=[sort
##' features or not], dataDist=[represents the present values],
##' inverse=[net is inverse mode or not], fCount=[amount of times
##' features were sorted], componentLayer=[which layer is the
##' 'bottleneck' (principal components)],
##' erro=[the used error function], gradient=[the used gradient method],
##' weights=[the present weights],
##' maxIter=[the amount of iterations that was done], scalingFactor=[the
##' scale of the original matrix])}
##' 
##' Slots
##' 
##' \describe{
##' \item{net}{"matrix",  matrix showing the representation of the
##' neural network, e.g. (2,4,6) for a network with two features, a
##' hidden layer and six output neurons (original variables).}
##' \item{hierarchic}{"list",  the hierarchic design of the network,
##' holds 'idx' (), 'var' () and layer (which layer is the principal
##' component layer).}
##' \item{fct}{"character",  a vector naming the functions that will be
##' applied on each layer. "linr" is linear (i.e.) standard matrix
##' products and "tanh" means that the arcus tangens is applied on the
##' result of the matrix product (for non-linearity).}
##' \item{fkt}{"character",  same as fct but the functions used during
##' back propagation.}
##' \item{weightDecay}{"numeric",  the value that is used to
##' incrementally decrease the weight changes to ensure convergence.}
##' \item{featureSorting}{"logical", indicates if features will be
##' sorted or not. This is used to make the NLPCA assume properties
##' closer to those of standard PCA were the first component is more
##' important for reconstructing the data than the second component.}
##' \item{dataDist}{"matrix", a matrix of ones and zeroes indicating
##' which values will add to the errror.}
##' \item{inverse}{"logical", network is inverse mode (currently only
##' inverse is supported) or not. Eg. the case when we have truly
##' missing values and wish to impute them.}
##' \item{fCount}{"integer", Counter for the amount of times features
##' were really sorted.}
##' \item{componentLayer}{"numeric", the index of 'net' that is the
##' component layer.}
##' \item{error}{"function", the used error function. Currently only one
##' is provided \code{errorHierarchic}.}
##' \item{gradient}{"function", the used gradient function. Currently
##' only one is provided  \code{derrorHierarchic}}
##' \item{weights}{"list", A list holding managements of the
##' weights. The list has two functions, weights$current() and
##' weights$set() which access a matrix in the local environment of
##' this object.}
##' \item{maxIter}{"integer", the amount of iterations used to train
##' this network.}
##' \item{scalingFactor}{"numeric", training the network is best made
##' with 'small' values so the original data is scaled down to a
##' suitable range by division with this number.}}
##' 
##' Methods
##'
##' \describe{ \item{vector2matrices}{Returns the
##' weights in a matrix representation.} }   
##' @title Class representation of the NLPCA neural net
##' @docType class
##' @aliases nlpcaNet nlpcaNet-class
##' @seealso \code{\link{nlpca}}
##' @aliases nFit nFit-class
##' @exportClass nlpcaNet
##' @keywords classes
##' @name pcaNet
##' @author Henning Redestig 
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
                   hierarchic=list(var=rbind(c(1,1,0.01)),
                     layer=3, idx=rbind(c(1,1,0),c(0,1,1))),
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
setAs("NULL", "nlpcaNet",
      function(from, to){
        new(to)
      })

##' This is a class representation of a PCA result
##' 
##' \bold{Creating Objects}\cr
##' \code{new("pcaRes", scores=[the scores], loadings=[the loadings],
##' nPcs=[amount of PCs], R2cum=[cumulative R2], nObs=[amount of
##' observations], nVar=[amount of variables], R2=[R2 for each
##' individual PC], sDev=[stdev for each individual PC],
##' centered=[was data centered], center=[original means],
##' varLimit=[what variance limit was exceeded], method=[method used to
##' calculate PCA], missing=[amount of NAs], 
##' completeObs=[estimated complete observations])}
##' 
##' \bold{Slots}\cr
##' \describe{
##'     \item{scores}{"matrix",  the calculated scores}
##'     \item{loadings}{"matrix",  the calculated loadings}
##'     \item{R2cum}{"numeric",  the cumulative R2 values}
##'     \item{sDev}{"numeric",  the individual standard
##'       deviations of the score vectors}
##'     \item{R2}{"numeric",  the individual R2 values}
##'     \item{cvstat}{"numeric",  cross-validation statistics}
##'     \item{nObs}{"numeric", number of observations}
##'     \item{nVar}{"numeric", number of variables}
##'     \item{centered}{"logical", data was centered or not}
##'     \item{center}{"numeric", the original variable centers}
##'     \item{scaled}{"logical", data was scaled or not}
##'     \item{scl}{"numeric", the original variable scales}
##'     \item{varLimit}{"numeric", the exceeded variance limit}
##'     \item{nPcs,nP}{"numeric", the number of calculated PCs}
##'     \item{method}{"character", the method used to perform PCA}
##'     \item{missing}{"numeric", the total amount of missing values in
##'       original data}
##'     \item{completeObs}{"matrix", the estimated complete observations}
##'     \item{network}{"nlpcaNet", the network used by non-linear PCA}
##'   }
##' 
##' \bold{Methods (not necessarily exhaustive)}\cr
##' \describe{
##'     \item{print}{Print function}
##'     \item{summary}{Extract information about PC relevance}
##'     \item{screeplot}{Plot a barplot of standard deviations for PCs}
##'     \item{slplot}{Make a side by side score and loadings plot}
##'     \item{nPcs}{Get the number of PCs}
##'     \item{nObs}{Get the number of observations}
##'     \item{cvstat}{Cross-validation statistics}
##'     \item{nVar}{Get the number of variables}
##'     \item{loadings}{Get the loadings}
##'     \item{scores}{Get the scores}
##'     \item{dim}{Get the dimensions (number of observations, number of
##'       features)}
##'     \item{centered}{Get a logical indicating if centering was done as
##'       part of the model}
##'     \item{center}{Get the averages of the original variables.}
##'     \item{completeObs}{Get the imputed data set}
##'     \item{method}{Get a string naming the used PCA method}
##'     \item{sDev}{Get the standard deviations of the PCs}
##'     \item{scaled}{Get a logical indicating if scaling was done as
##'       part of the model}
##'     \item{scl}{Get the scales of the original variablesb}
##'     \item{R2cum}{Get the cumulative R2}
##'   }
##' @title Class for representing a PCA result
##' @keywords classes
##' @exportClass pcaRes
##' @docType class
##' @name pcaRes
##' @aliases pcaRes pcaRes-class 
##' @author Henning Redestig
setClass("pcaRes",
         representation(completeObs="matrix",
			scores="matrix",
                        loadings="matrix",
                        R2cum="numeric",
                        R2="numeric",   # ditch, get from R2cum
                        cvstat="numeric",   # ditch, get from R2cum
                        sDev="numeric", # ditch, get from scores
                        nObs="numeric", # ditch, get from scores
                        nVar="numeric", 
                        centered="logical", 
                        center="numeric",
                        subset="numeric",
                        scaled="character", 
                        scale="numeric",
                        varLimit="numeric", # ditch, useless
                        nPcs="numeric",     # ditch, get from scores
                        method="character",
                        missing="matrix",
                        network="nlpcaNet"),
         prototype(completeObs=NULL,
		   scores=NULL,
                   loadings=NULL,
                   R2cum=NULL,
                   R2=NULL,
                   subset=NULL,
                   cvstat=NULL,
                   sDev=NULL,
                   nObs=NULL,
                   nVar=NULL,
                   centered=NULL,
                   center=NULL,
                   scaled=NULL,
                   scale=NULL,
                   varLimit=NULL,
                   nPcs=NULL,
                   method=NULL,
                   missing=NULL,
                   network=NULL))
setAs("NULL", "pcaRes",
      function(from, to){
        new(to)
      })

##' This is a class representation of nearest neighbour imputation
##' (nni) result
##'
##' \bold{Creating Objects}\cr
##' \code{new("nniRes", completeObs=[the estimated complete
##' observations], k=[cluster size], nObs=[amount of observations],
##' nVar=[amount of variables], centered=[was the data centered befor
##' running LLSimpute], center=[original means], method=[method used
##' to perform clustering], missing=[amount of NAs])}
##' 
##' \bold{Slots}\cr
##' \describe{
##'   \item{completeObs}{"matrix", the estimated complete observations}
##'   \item{nObs}{"numeric", amount of observations}
##'   \item{nVar}{"numeric", amount of variables}
##'   \item{correlation}{"character", the correlation method used
##'         (pearson, kendall or spearman)}
##'   \item{centered}{"logical", data was centered or not}
##'   \item{center}{"numeric", the original variable centers}
##'   \item{k}{"numeric", cluster size}
##'   \item{method}{"character", the method used to perform the clustering}
##'   \item{missing}{"numeric", the total amount of missing values in
##'     original data}
##' }
##' 
##' \bold{Methods}\cr
##' \describe{ \item{print}{Print function} }
##' @title Class for representing a nearest neighbour imputation result
##' @docType class
##' @exportClass nniRes
##' @name nniRes
##' @keywords classes
##' @aliases nniRes nniRes-class
##' @author Wolfram Stacklies
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

##' Create an object that holds the weights for nlpcaNet. Holds and
##' sets weights in using an environment object.
##' @param w  \code{matrix} -- New weights
##' @return A weightsAccound with \code{set} and \code{current}
##' functions.
##' @author Henning Redestig
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

