##' Simulated data set looking like a helix
##' 
##' 
##' A matrix containing 1000 observations (rows) and three variables
##' (columns).
##' @title A helix structured toy data set
##' @name helix
##' @aliases helix
##' @usage data(helix)
##' @docType data
##' @references Matthias Scholz, Fatma Kaplan, Charles L. Guy, Joachim
##' Kopka and Joachim Selbig. - Non-linear PCA: a missing data
##' approach.  \emph{Bioinformatics 2005 21(20):3887-3895}
##' @keywords datasets
##' @author Henning Redestig 
NULL

##' A complete subset from a larger metabolite data set.  This is the
##' original, complete data set and can be used to compare estimation
##' results created with the also provided incomplete data (called
##' metaboliteData).  The data was created during an in house
##' Arabidopsis coldstress experiment.
##'
##' A matrix containing 154 observations (rows) and 52 metabolites
##' (columns).
##' @name metaboliteDataComplete
##' @docType data
##' @aliases metaboliteDataComplete
##' @title A complete metabolite data set from an Arabidopsis
##' coldstress experiment
##' @keywords datasets
##' @seealso  \code{\link{metaboliteData}}
##' @references Matthias Scholz, Fatma Kaplan, Charles L. Guy, Joachim
##' Kopka and Joachim Selbig. - Non-linear PCA: a missing data
##' approach.\emph{Bioinformatics 2005 21(20):3887-3895}
##' @author Wolfram Stacklies
NULL

##' A incomplete subset from a larger metabolite data set.  This is the
##' original, complete data set and can be used to compare estimation
##' results created with the also provided incomplete data (called
##' metaboliteData).  
##'
##' A matrix containing 154 observations (rows) and 52 metabolites
##' (columns). The data contains 5\% of artificially created uniformly
##' distributed misssing values. The data was created during an in
##' house Arabidopsis coldstress experiment.
##' @name metaboliteData
##' @docType data
##' @aliases metaboliteData
##' @title A incomplete metabolite data set from an Arabidopsis
##' coldstress experiment
##' @keywords datasets
##' @seealso  \code{\link{metaboliteDataComplete}}
##' @references Matthias Scholz, Fatma Kaplan, Charles L. Guy, Joachim
##' Kopka and Joachim Selbig. - Non-linear PCA: a missing data
##' approach.\emph{Bioinformatics 2005 21(20):3887-3895}
##' @author Wolfram Stacklies
NULL


##' Principal Component Analysis in R
##'
##' \tabular{ll}{
##' Package: \tab pcaMethods \cr
##' Type: \tab Package \cr
##' Developed since: \tab 2006 \cr
##' License: \tab GPL (>=3) \cr
##' LazyLoad: \tab yes \cr
##' }
##'
##' Provides Bayesian PCA, Probabilistic PCA, Nipals PCA, Inverse
##' Non-Linear PCA and the conventional SVD PCA. A cluster  based
##' method for missing value estimation is included for comparison.
##' BPCA, PPCA and NipalsPCA may be used to perform PCA on incomplete
##' data as well as for accurate missing value estimation.  A set of
##' methods for printing and plotting the results is also provided.
##' All PCA methods make use of the same data structure (pcaRes) to
##' provide a unique interface to the PCA results. Developed at the
##' Max-Planck Institute for Molecular Plant Physiology, Golm,
##' Germany, RIKEN Plant Science Center Yokohama, Japan, and CAS-MPG
##' Partner Institute for Computational Biology (PICB) Shanghai,
##' P.R. China
##'
##' @name pcaMethods
##' @aliases pcaMethods
##' @docType package
##' @importFrom Rcpp evalCpp
##' @import Biobase
##' @import BiocGenerics
##' @import methods
##' @title pcaMethods
##' @useDynLib pcaMethods
##' @author Wolfram Stacklies, Henning Redestig
NULL

##' \describe{
##' \item{plotR2}{Lack of relevance for this plot and the fact that it
##' can not show cross-validation based diagnostics in the same plot
##' makes it redundant with the introduction of a dedicated
##' \code{plot} function for  \code{pcaRes}. The new plot only shows
##' R2cum but the result is pretty much the same.}}
##' @name pcaMethods-deprecated
##' @aliases pcaMethods-deprecated
##' @title Deprecated methods for pcaMethods
##' @author Henning Redestig 
NULL
