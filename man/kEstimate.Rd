\encoding{latin1}
\name{KEstimate}
\alias{kEstimate}
\title{Estimate best number of Components for missing value estimation}
\description{
  Perform cross validation to estimate the optimal number of components for missing
  value estimation.
  Cross validation is done for the complete subset of a variable.
  The assumption hereby is that variables that are highly correlated in a distinct region (here
  the non-missing observations) are also correlated in another (here the missing observations).
  This also implies that the complete subset must be large enough to be representative.
  For each incomplete variable, the available values are divided into a user defined
  number of cv-segments. The segments have equal size, but are chosen from a random
  equal distribution. The non-missing values of the variable are covered completely.
  PPCA, BPCA, SVDimpute, Nipals PCA, llsImpute an NLPCA may be used for imputation. \cr
  The whole cross validation is repeated several times so, depending on
  the parameters, the calculations can take very long time.
  As error measure the NRMSEP (see Feten et. al, 2005) or the Q2 distance is used.
  The NRMSEP basically
  normalises the RMSD between original data and estimate by the variable-wise
  variance. The reason for this is that a higher variance will generally lead to a
  higher estimation error.
  If the number of samples is small, the gene - wise variance may become an unstable
  criterion and the Q2 distance should be used instead. Also if variance normalisation
  was applied previously.
}
\details{
  Run time may be very high on large data sets. Especially when used with methods
  complex methods like BPCA or Nipals PCA.
  The estimation method is called 
  \eqn{(g_{miss} \cdot segs \cdot nruncv)}{(g\_miss * segs * nruncv)}
  times, where \eqn{g_{miss}}{g\_miss}
  is the number of variables showing missing values.\cr
  As cross validation can only be done gene-wise, in this function Q2 is
  defined on single genes, not on the entire data set. This is when calculating Q2
  as \eqn{\frac{\sum(X - Xe)^2}{\sum(X^2)}}{sum(X - Xe)^2 \ sum(X^2)},
  X is the currently used gene and Xe it's estimate.
}
\usage{
kEstimate(Matrix, method = "ppca", evalPcs = 1:3, segs = 3, nruncv = 10,
em = "nrmsep", allVariables = FALSE, verbose = interactive(), 
random = FALSE,...)
}
\arguments{
\item{Matrix}{\code{matrix} -- numeric matrix containing observations in rows and 
  variables in columns}
\item{method}{\code{character} -- One of ppca | bpca | svdImpute | nipals | nlpca | llsImpute
  | llsImputeAll. llsImputeAll uses llsImpute with option allVariables = TRUE.}
\item{evalPcs}{\code{numeric} -- The principal components to use for cross validation
or cluster sizes if used with llsImpute.
Should be an array containing integer values, eg. evalPcs = 1:10
or evalPcs = C(2,5,8).The NRMSEP is calculated for each component.}
\item{segs}{\code{numeric} -- number of segments for cross validation}
\item{nruncv}{\code{numeric} -- Times the whole cross validation is repeated}
\item{em}{\code{character} -- The error measure. This can be nrmsep or q2}
\item{allVariables}{\code{boolean} -- If TRUE, the NRMSEP is calculated for all variables,
If FALSE, only the incomplete ones are included.
You maybe want to do this to compare several methods on a 
complete data set.}
\item{verbose}{\code{boolean} -- If TRUE, the NRMSEP and the variance are printed
to the console each iteration.}
\item{random}{\code{boolean} -- Impute normal distributed random values with
same mean and standard deviation than the original data.
This is only thought for comparison.}
\item{...}{Further arguments to \code{pca}}
}
\value{
  \item{list}{Returns a list with the elements:
  \itemize{
    \item mink - number of PCs for which the minimal average NRMSEP was obtained
    \item eError - an array of of size length(evalPcs). Contains the average error
    of the cross validation runs for each number of components.
    \item eErrorMat - Estimation error for each of the separately done cross validations.
    Rows represent the separate runs, columns the number of components
    used for estimation.
    \item evalPcs - The evaluated numbers of components or cluster sizes 
    (the same as the evalPcs input parameter).
  }}
}
\seealso{
  \code{\link{bpca}, \link{svdImpute}, \link{prcomp}, \link{nipalsPca}, \link{pca}}.
}
\examples{
## Load a sample metabolite dataset (metaboliteData)
data(metaboliteData)

# Now remove 10\% of the data
cond<-runif(length(metaboliteData))
cond<-cond<0.1
metaboliteData[cond] <- NA

# Do cross validation with ppca for component 2:4
esti <- kEstimate(metaboliteData, method = "ppca", evalPcs = 2:4, nruncv=1, em="nrmsep")

# Plot the result
barplot(drop(esti$eError), xlab = "Components",ylab = "NRMSEP (1 iterations)")

# The best k value is:
print(esti$mink)
}
\keyword{multivariate}

\author{Wolfram Stacklies \cr
        CAS-MPG Partner Institute for Computational Biology, Shanghai, China \cr        
        \email{wolfram.stacklies@gmail.com} \cr
}

