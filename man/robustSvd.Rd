\name{robustSvd}
\alias{robustSvd}
\title{Alternating L1 Singular Value Decomposition}
\usage{
  robustSvd(x)
}
\arguments{
  \item{x}{A matrix whose SVD decomposition is to be
  computed. Missing values are allowed.}
}
\value{
  The robust SVD of the matrix is x = u d v'. \item{d}{A
  vector containing the singular values of \code{x}.}
  \item{u}{A matrix whose columns are the left singular
  vectors of \code{x}.} \item{v}{A matrix whose columns are
  the right singular vectors of \code{x}.}
}
\description{
  A robust approximation to the singular value
  decomposition of a rectangular matrix is computed using
  an alternating L1 norm (instead of the more usual least
  squares L2 norm). As the SVD is a least-squares
  procedure, it is highly susceptible to outliers and in
  the extreme case, an individual cell (if sufficiently
  outlying) can draw even the leading principal component
  toward itself.
}
\details{
  See Hawkins et al (2001) for details on the robust SVD
  algorithm. Briefly, the idea is to sequentially estimate
  the left and right eigenvectors using an L1 (absolute
  value) norm minimization.

  Note that the robust SVD is able to accomodate missing
  values in the matrix \code{x}, unlike the usual
  \code{svd} function.

  Also note that the eigenvectors returned by the robust
  SVD algorithm are NOT (in general) orthogonal and the
  eigenvalues need not be descending in order.
}
\note{
  Two differences from the usual SVD may be noted. One
  relates to orthogonality. In the conventional SVD, all
  the eigenvectors are orthogonal even if not explicitly
  imposed.  Those returned by the AL1 algorithm (used here)
  are (in general) not orthogonal. Another difference is
  that, in the L2 analysis of the conventional SVD, the
  successive eigen triples (eigenvalue, left eigenvector,
  right eigenvector) are found in descending order of
  eigenvalue. This is not necessarily the case with the AL1
  algorithm.  Hawkins et al (2001) note that a larger eigen
  value may follow a smaller one.
}
\examples{
## Load a complete sample metabolite data set and mean center the data
data(metaboliteDataComplete)
mdc <- prep(metaboliteDataComplete, center=TRUE, scale="none")
## Now create 5\% of outliers.
cond   <- runif(length(mdc)) < 0.05;
mdcOut <- mdc
mdcOut[cond] <- 10
## Now we do a conventional SVD and a robustSvd on both, the original and the
## data with outliers.
resSvd <- svd(mdc)
resSvdOut <- svd(mdcOut)
resRobSvd <- robustSvd(mdc)
resRobSvdOut <- robustSvd(mdcOut)
## Now we plot the results for the original data against those with outliers
## We can see that robustSvd is hardly affected by the outliers.
plot(resSvd$v[,1], resSvdOut$v[,1])
plot(resRobSvd$v[,1], resRobSvdOut$v[,1])
}
\author{
  Kevin Wright, modifications by Wolfram Stacklies
}
\references{
  Hawkins, Douglas M, Li Liu, and S Stanley Young (2001)
  Robust Singular Value Decomposition, National Institute
  of Statistical Sciences, Technical Report Number 122.
  \url{http://www.niss.org/technicalreports/tr122.pdf}
}
\seealso{
  \code{\link{svd}}, \code{\link[ade4:nipals]{nipals}} for
  an alternating L2 norm method that also accommodates
  missing data.
}
\keyword{algebra}

