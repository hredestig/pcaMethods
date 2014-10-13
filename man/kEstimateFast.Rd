\name{kEstimateFast}
\alias{kEstimateFast}
\title{Estimate best number of Components for missing value estimation}
\usage{
  kEstimateFast(Matrix, method = "ppca", evalPcs = 1:3,
    em = "nrmsep", allVariables = FALSE,
    verbose = interactive(), ...)
}
\arguments{
  \item{Matrix}{\code{matrix} -- numeric matrix containing
  observations in rows and variables in columns}

  \item{method}{\code{character} -- a valid pca method (see
  \code{\link{pca}}).}

  \item{evalPcs}{\code{numeric} -- The principal components
  to use for cross validation or cluster sizes if used with
  llsImpute. Should be an array containing integer values,
  eg. evalPcs = 1:10 or evalPcs = C(2,5,8).The NRMSEP is
  calculated for each component.}

  \item{em}{\code{character} -- The error measure. This can
  be nrmsep or q2}

  \item{allVariables}{\code{boolean} -- If TRUE, the NRMSEP
  is calculated for all variables, If FALSE, only the
  incomplete ones are included. You maybe want to do this
  to compare several methods on a complete data set.}

  \item{verbose}{\code{boolean} -- If TRUE, the NRMSEP and
  the variance are printed to the console each iteration.}

  \item{...}{Further arguments to \code{pca}}
}
\value{
  \item{list}{Returns a list with the elements: \itemize{
  \item minNPcs - number of PCs for which the minimal
  average NRMSEP was obtained \item eError - an array of of
  size length(evalPcs). Contains the estimation error for
  each number of components. \item evalPcs - The evaluated
  numbers of components or cluster sizes (the same as the
  evalPcs input parameter). }}
}
\description{
  This is a simple estimator for the optimal number of
  componets when applying PCA or LLSimpute for missing
  value estimation.  No cross validation is performed,
  instead the estimation quality is defined as
  Matrix[!missing] - Estimate[!missing]. This will give a
  relatively rough estimate, but the number of iterations
  equals the length of the parameter evalPcs.\cr Does not
  work with LLSimpute!! As error measure the NRMSEP (see
  Feten et. al, 2005) or the Q2 distance is used.  The
  NRMSEP basically normalises the RMSD between original
  data and estimate by the variable-wise variance. The
  reason for this is that a higher variance will generally
  lead to a higher estimation error.  If the number of
  samples is small, the gene - wise variance may become an
  unstable criterion and the Q2 distance should be used
  instead. Also if variance normalisation was applied
  previously.
}
\examples{
data(metaboliteData)
# Estimate best number of PCs with ppca for component 2:4
esti <- kEstimateFast(t(metaboliteData), method = "ppca", evalPcs = 2:4, em="nrmsep")
barplot(drop(esti$eError), xlab = "Components",ylab = "NRMSEP (1 iterations)")
# The best k value is:
print(esti$minNPcs)
}
\author{
  Wolfram Stacklies
}
\seealso{
  \code{\link{kEstimate}}.
}
\keyword{multivariate}

