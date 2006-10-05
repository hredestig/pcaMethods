\encoding{latin1}
\name{svdImpute}
\alias{svdImpute}
\title{SVDimpute algorithm}
\description{
This implements the SVDimpute algorithm as proposed by Troyanskaya et al, 2001.
The idea behind the algorithm is to estimate the missing values as a
linear combination of the \code{k} most significant eigengenes.

Missing values are denoted as \code{NA}
}
\details{
  As SVD can only be performed on complete matrices, all missing
  values are initially replaced by 0 (what is in fact the mean
  on centred data).
  The algorithm works iteratively until the change in the estimated solution
  falls below a certain threshold.
  Each step the eigengenes of the current estimate are calculated and
  used to determine a new estimate. Eigengenes denote the loadings if
  pca is performed considering genes as observations.

  An optimal linear combination is found by regressing the incomplete
  gene against the \code{k} most significant eigengenes. If the value at
  position \code{j} is missing,
  the \eqn{j^th}{j^th} value of the eigengenes is not used when determining
  the regression coefficients.

  \bold{Complexity:} Each iteration, standard PCA (\code{prcomp}) needs to be
  done for each incomplete gene to get the eigengenes. This is usually fast for
  small data sets, but complexity may rise if the data sets become very large.
}
\usage{
svdImpute(Matrix, nPcs = 2, center=TRUE, completeObs=TRUE, threshold = 0.01, 
  maxSteps = 100, verbose = interactive(), ...)
}
\arguments{
  \item{Matrix}{\code{matrix} -- Data containing the variables in
    columns and observations in rows. The data may contain missing values,
    denoted as \code{NA}.}
  \item{nPcs}{\code{numeric} -- Number of components to estimate.
    The preciseness of the missing value estimation depends on the
    number of components, which should resemble the internal structure
    of the data.}
  \item{center}{Mean center the data if TRUE}
  \item{completeObs}{Return the estimated complete observations if TRUE. This is
    the input data with NA values replaced by the estimated values.}
  \item{threshold}{The iteration stops if the change in the matrix
    falls below this threshold, the default is 0.01.
    (0.01 was empirically determined by Troyanskaya et. al)}
  \item{maxSteps}{Maximum number of iteration steps. Default is 100.}
  \item{verbose}{Print some output if TRUE. Default is interactive()}
  \item{...}{Reserved for parameters used in future version of the algorithm}
}
\value{
  \item{pcaRes}{Standart PCA result object used by all
    PCA-based methods of this package. Contains scores, loadings, data mean and
    more. See \code{\link{pcaRes}} for details.}
}
\seealso{
  \code{\link{bpca}, \link{ppca}, \link{prcomp}, \link{nipalsPca}, \link{pca}, \link{pcaRes}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData)
data(metaboliteData)

# Now remove 10\% of the data
rows <- nrow(metaboliteData)
cols <- ncol(metaboliteData)
cond<-matrix(runif(rows * cols),rows,cols) < 0.1
metaboliteData[cond] <- NA

## Perform probabilistic PCA using the 3 largest components
result <- pca(metaboliteData, method="svdImpute", nPcs=3, center = TRUE)

## Get the estimated principal axes (loadings)
loadings <- result@loadings

## Get the estimated scores
scores <- result@scores

## Get the estimated complete observations
cObs <- result@completeObs

## Now plot the scores
plotPcs(result, scoresLoadings=c(TRUE,FALSE))

}
\keyword{multivariate}

\author{Wolfram Stacklies \cr
        Max Planck Institut fuer Molekulare Pflanzenphysiologie, Potsdam, Germany \cr
        \email{wolfram.stacklies@gmail.com} \cr
}
\references{
Troyanskaya O. and Cantor M. and Sherlock G. and Brown
P. and Hastie T. and Tibshirani R. and Botstein D. and Altman RB.
- Missing value estimation methods for DNA microarrays.
\emph{Bioinformatics. 2001 Jun;17(6):520-5.}
}
