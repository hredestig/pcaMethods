\name{robustPca}
\alias{robustPca}
\title{PCA implementation based on robustSvd}
\usage{
  robustPca(Matrix, nPcs = 2, verbose = interactive(), ...)
}
\arguments{
  \item{Matrix}{\code{matrix} -- Data containing the
  variables in columns and observations in rows. The data
  may contain missing values, denoted as \code{NA}.}

  \item{nPcs}{\code{numeric} -- Number of components to
  estimate. The preciseness of the missing value estimation
  depends on the number of components, which should
  resemble the internal structure of the data.}

  \item{verbose}{\code{boolean} Print some output to the
  command line if TRUE}

  \item{...}{Reserved for future use. Currently no further
  parameters are used}
}
\value{
  Standard PCA result object used by all PCA-based methods
  of this package. Contains scores, loadings, data mean and
  more. See \code{\link{pcaRes}} for details.  are used.
}
\description{
  This is a PCA implementation robust to outliers in a data
  set. It can also handle missing values, it is however NOT
  intended to be used for missing value estimation.  As it
  is based on robustSVD we will get an accurate estimation
  for the loadings also for incomplete data or for data
  with outliers.  The returned scores are, however,
  affected by the outliers as they are calculated inputData
  X loadings. This also implies that you should look at the
  returned R2/R2cum values with caution.  If the data show
  missing values, scores are caluclated by just setting all
  NA - values to zero. This is not expected to produce
  accurate results. Please have also a look at the manual
  page for \code{robustSvd}. Thus this method should mainly
  be seen as an attempt to integrate \code{robustSvd()}
  into the framework of this package.  Use one of the other
  methods coming with this package (like PPCA or BPCA) if
  you want to do missing value estimation.  It is not
  recommended to use this function directely but rather to
  use the pca() wrapper function.
}
\details{
  The method is very similar to the standard
  \code{prcomp()} function.  The main difference is that
  \code{robustSvd()} is used instead of the conventional
  \code{svd()} method.
}
\examples{
## Load a complete sample metabolite data set and mean center the data
data(metaboliteDataComplete)
mdc <- scale(metaboliteDataComplete, center=TRUE, scale=FALSE)
## Now create 5\\\% of outliers.
cond   <- runif(length(mdc)) < 0.05;
mdcOut <- mdc
mdcOut[cond] <- 10
## Now we do a conventional PCA and robustPca on the original and the data
## with outliers.
## We use center=FALSE here because the large artificial outliers would
## affect the means and not allow to objectively compare the results.
resSvd    <- pca(mdc, method = "svd", nPcs = 10, center = FALSE)
resSvdOut <- pca(mdcOut, method = "svd", nPcs = 10, center = FALSE)
resRobPca <- pca(mdcOut, method = "robustPca", nPcs = 10, center = FALSE)
## Now we plot the results for the original data against those with outliers
## We can see that robustPca is hardly effected by the outliers.
plot(loadings(resSvd)[,1], loadings(resSvdOut)[,1])
plot(loadings(resSvd)[,1], loadings(resRobPca)[,1])
}
\author{
  Wolfram Stacklies
}
\seealso{
  \code{\link{robustSvd}, \link{svd}, \link{prcomp},
  \link{pcaRes}}.
}
\keyword{multivariate}

