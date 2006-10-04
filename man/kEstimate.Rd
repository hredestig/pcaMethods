\encoding{latin1}
\name{KEstimate}
\alias{kEstimate}
\title{Estimate best number of Components for missing value estimation}
\description{
  Perform cross validation to estimate the optimal number of components for missing
  value estimation.
  Cross validation is done on the subset containing only complete observations
  because including incomplete observations may tamper the results.
  The assumption hereby is that genes that are highly correlated in a distinct region (here
  the non-missing observations) are also correlated in another (here the missing observations).
  This also implies that the complete subset must be large enough to be representative.
  For each incomplete gene, the available values are divided into a user defined
  number of cv-segments. The segments have equal size, but are chosen from a random
  equal distribution. The non-missing values of the gene are covered completely.
  PPCA, BPCA, SVDimpute and Nipals PCA may be used for imputation. \cr
  The whole cross validation is repeated several times.
  As error measure the NRMSEP (see Feten et. al, 2005) is used. This error
  basically
  normalises the RMSD between original data and estimate by the gene-wise
  variance. The reason for this is that a higher variance will lead to a
  higher estimation error.
}
\details{
  Run time may be very high on large data sets. Also, when used with methods
  like BPCA or Nipals PCA which are already quite slow.
  The estimation method is called 
  \eqn{(g_{miss} \cdot segs \cdot nruncv)}{(g\_miss * segs * nruncv)}
  times, where \eqn{g_{miss}}{g\_miss}
  is the number of genes showing missing values.
}
\usage{
kEstimate(data, method = "ppca", maxPcs = 3, segs = 3, nruncv = 10,
allGenes = FALSE, verbose = interactive(), random = FALSE)
}
\arguments{
\item{data}{\code{matrix} -- numeric matrix containing observations in rows and 
  genes in columns}
\item{method}{\code{character} -- One of ppca | bpca | svdImpute | nipals}
\item{maxPcs}{\code{numeric} -- number of principal components to use for cross validation.
The NRMSEP is calculated for 1:maxPcs components.}
\item{segs}{\code{numeric} -- number of segments for cross validation}
\item{nruncv}{\code{numeric} -- Times the whole cross validation is repeated}
\item{allGenes}{\code{boolean} -- If TRUE, the NRMSEP is calculated for all genes,
If FALSE, only the incomplete ones are included.
You maybe want to do this to compare several methods on a 
complete data set.}
\item{verbose}{\code{boolean} -- If TRUE, the NRMSEP and the variance are printed
to the console each iteration.}
\item{random}{\code{boolean} -- Impute normal distributed random values with
same mean and standard deviation than the original data.
This is only thought for comparison.}
}
\value{
  \item{list}{Returns a list with the elements:
  \itemize{
    \item mink - number of PCs for which the minimal average NRMSEP was obtained
    \item nrmsep - a matrix of dimension (nruncv, maxPcs).
    The columns contain the NRMSEP obtained for each repeat of the cross validation.
  }}
}
\seealso{
  \code{\link{bpca}, \link{svdImpute}, \link{prcomp}, \link{nipalsPca}, \link{pca}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData)
data(metaboliteData)

# Now remove 10\% of the data
rows <- nrow(metaboliteData)
cols <- ncol(metaboliteData)
cond<-matrix(runif(rows * cols),rows,cols) < 0.1
metaboliteData[cond] <- NA

# Do cross validation with ppca for component 1:3
nrmsep <- kEstimate(metaboliteData, method = "ppca", maxPcs = 3, nruncv=1)

# Plot the result
barplot(drop(nrmsep$nrmsep), xlab = "Components",ylab = "NRMSEP (1 iterations)")

# The best k value is:
print(nrmsep$mink)
}
\keyword{multivariate}

\author{Wolfram Stacklies \cr
        Max Planck Institut fuer Molekulare Pflanzenphysiologie, Potsdam, Germany \cr
        \email{wolfram.stacklies@gmail.com} \cr
}

