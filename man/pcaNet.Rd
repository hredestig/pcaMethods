\docType{class}
\name{pcaNet}
\alias{nFit}
\alias{nFit-class}
\alias{nlpcaNet}
\alias{nlpcaNet-class}
\alias{pcaNet}
\title{Class representation of the NLPCA neural net}
\description{
  This is a class representation of a non-linear PCA neural
  network. The \code{nlpcaNet} class is not meant for
  user-level usage.
}
\details{
  Creating Objects

  \code{new("nlpcaNet", net=[the network structure],
  hierarchic=[hierarchic design], fct=[the functions at
  each layer], fkt=[the functions used for forward
  propagation], weightDecay=[incremental decrease of weight
  changes over iterations (between 0 and 1)],
  featureSorting=[sort features or not],
  dataDist=[represents the present values], inverse=[net is
  inverse mode or not], fCount=[amount of times features
  were sorted], componentLayer=[which layer is the
  'bottleneck' (principal components)], erro=[the used
  error function], gradient=[the used gradient method],
  weights=[the present weights], maxIter=[the amount of
  iterations that was done], scalingFactor=[the scale of
  the original matrix])}

  Slots

  \describe{ \item{net}{"matrix", matrix showing the
  representation of the neural network, e.g. (2,4,6) for a
  network with two features, a hidden layer and six output
  neurons (original variables).} \item{hierarchic}{"list",
  the hierarchic design of the network, holds 'idx' (),
  'var' () and layer (which layer is the principal
  component layer).} \item{fct}{"character", a vector
  naming the functions that will be applied on each layer.
  "linr" is linear (i.e.) standard matrix products and
  "tanh" means that the arcus tangens is applied on the
  result of the matrix product (for non-linearity).}
  \item{fkt}{"character", same as fct but the functions
  used during back propagation.}
  \item{weightDecay}{"numeric", the value that is used to
  incrementally decrease the weight changes to ensure
  convergence.} \item{featureSorting}{"logical", indicates
  if features will be sorted or not. This is used to make
  the NLPCA assume properties closer to those of standard
  PCA were the first component is more important for
  reconstructing the data than the second component.}
  \item{dataDist}{"matrix", a matrix of ones and zeroes
  indicating which values will add to the errror.}
  \item{inverse}{"logical", network is inverse mode
  (currently only inverse is supported) or not. Eg. the
  case when we have truly missing values and wish to impute
  them.} \item{fCount}{"integer", Counter for the amount of
  times features were really sorted.}
  \item{componentLayer}{"numeric", the index of 'net' that
  is the component layer.} \item{error}{"function", the
  used error function. Currently only one is provided
  \code{errorHierarchic}.} \item{gradient}{"function", the
  used gradient function. Currently only one is provided
  \code{derrorHierarchic}} \item{weights}{"list", A list
  holding managements of the weights. The list has two
  functions, weights$current() and weights$set() which
  access a matrix in the local environment of this object.}
  \item{maxIter}{"integer", the amount of iterations used
  to train this network.} \item{scalingFactor}{"numeric",
  training the network is best made with 'small' values so
  the original data is scaled down to a suitable range by
  division with this number.}}

  Methods

  \describe{ \item{vector2matrices}{Returns the weights in
  a matrix representation.} }
}
\author{
  Henning Redestig
}
\seealso{
  \code{\link{nlpca}}
}
\keyword{classes}

