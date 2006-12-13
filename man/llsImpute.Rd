\encoding{latin1}
\name{llsImpute}
\alias{llsImpute}
\title{LLSimpute algorithm}
\description{
  Missing value estimation using local least squares (LLS).
  First, k genes are selected by pearson, spearman or kendall correlation coefficients.
  Then missing values are imputed by a linear combination of the k selected
  genes. The optimal combination is found by LLS regression.
  The method was first described by Kim et al, Bioinformatics, 21(2),2005.

  Missing values are denoted as \code{NA}

  It is not recommended to use this function directely but rather to use the
  nni() wrapper function.
}
\details{
  The methods provides two ways for missing value estimation, selected by
  the \code{allGenes} option. The first one is to use only complete genes for the 
  regression. This is preferable when the number of incomplete genes is relatively small.

  The second way is to consider all genes as candidates for the regression.
  Hereby missing values are initially replaced by the columns wise mean. 
  The method then iterates, using
  the current estimate as input for the regression until the change between new
  and old estimate falls below a threshold (0.001).

  \bold{Complexity:} Each step the generalized inverse of a \code{miss} x {k}
  matrix is calculated. Where \code{miss} is the number of missing values in
  gene j and \code{k} the number of neighbours. This may be slow for large values of
  k and / or many missing values. See also help("ginv").
}
\usage{
  llsImpute(Matrix, k = 10, center = FALSE, completeObs = TRUE, correlation = "pearson",
  allGenes = FALSE, maxSteps = 100, verbose = interactive(), ...)
}
\arguments{
  \item{Matrix}{\code{matrix} -- Data containing the variables (genes) in
    columns and observations (samples) in rows. The data may contain missing values,
    denoted as \code{NA}.}
  \item{k}{\code{numeric} -- Cluster size, this is the number of similar genes
    used for regression.}
  \item{center}{\code{boolean} -- Mean center the data if TRUE}
  \item{completeObs}{\code{boolean} -- Return the estimated complete observations if 
    TRUE. This is the input data with NA values replaced by the estimated values.}
  \item{correlation}{\code{character} -- How to calculate the distance between genes. 
     One out of pearson | kendall | spearman , see also help("cor").}
  \item{allGenes}{\code{boolean} -- Use only complete genes to do the regression if TRUE, all
     genes if FALSE.}
  \item{maxSteps}{\code{numeric} -- Maximum number of iteration steps if allGenes = TRUE.}
  \item{verbose}{\code{boolean} -- Print step number and relative change if TRUE and 
    allGenes = TRUE}
  \item{...}{Reserved for parameters used in future version of the algorithm}
}
\value{
  \item{nniRes}{Standard nni (nearest neighbour imputation) result object of this package.
    See \code{\link{nniRes}} for details.}
}
\seealso{
  \code{\link{pca}, \link{nniRes}, \link{nni}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData) with already 10\% of
## data missing
data(metaboliteData)

## Perform llsImpute using k = 10
result <- llsImpute(metaboliteData, k = 10, correlation = "pearson")

## Get the estimated complete observations
cObs <- result@completeObs

}
\keyword{multivariate}

\author{Wolfram Stacklies \cr
        MPG/CAS Partner Institute for Computational Biology, Shanghai, P.R. China \cr
        \email{wolfram.stacklies@gmail.com} \cr
}
\references{
Kim, H. and Golub, G.H. and Park, H.
- Missing value estimation for DNA microarray gene expression data:
local least squares imputation.
\emph{Bioinformatics, 2005; 21(2):187-198.}

Troyanskaya O. and Cantor M. and Sherlock G. and Brown
P. and Hastie T. and Tibshirani R. and Botstein D. and Altman RB.
- Missing value estimation methods for DNA microarrays.
\emph{Bioinformatics. 2001 Jun;17(6):520-525.}
}
