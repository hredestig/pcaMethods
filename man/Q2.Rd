\name{Q2}
\alias{Q2}
\title{Cross-validation for PCA}
\usage{Q2(object, originalData, fold=5, nruncv=1, type=c("krzanowski",
    "impute"), verbose=interactive(), ...)}
\description{Internal cross-validation can be used for estimating the
level of structure in a data set and to optimise the choice of number
of principal components.}
\details{This method calculates \eqn{Q^2} for a PCA model. This is the predictory
version of \eqn{R^2} and can be interpreted as the ratio of
variance that can be predicted independetly by the PCA  model. Poor (low) \eqn{Q^2}
indicates that the PCA model only describes noise
and that the model is unrelated to the true data structure. The
definition of \eqn{Q^2} is:
\deqn{Q^2 = 1 - \frac{\sum_{i}^{k}\sum_{j}^{n}(x -
\hat{x})^2}{\sum_{i}^{k}\sum_{j}^{n}x^2}}{Q^2 = 1 - sum_i^k sum_j^n (x -
\hat{x})^2 / \sum_i^k \sum_j^n(x^2)}
for the matrix \eqn{x} which has \eqn{n} rows and \eqn{k}
columns. For a given amount of PC's x is estimated as \eqn{\hat{x} =
TP'} (T are scores and P are loadings). Although this defines the
leave-one-out cross-validation this is 
not what is performed if fold is less than the amount of rows and/or
columns.
In 'impute' type CV, diagonal rows of elements in the matrix are deleted and the
re-estimated. 
In 'krzanowski' type CV, rows are sequentially left out to build fold PCA models
which give the loadings. Then, columns are sequentially left out to build
fold models for scores. By combining scores and loadings from different models,
we can estimate completely left out values.
The two types may seem similar but can give very different results,
krzanowski typically yields more stable and reliable result for estimating data structure
whereas impute is better for evaluating missing value imputation performance.}
\value{A matrix or vector with \eqn{Q^2} estimates.}
\author{Henning Redestig}
\keyword{multivariate}
\arguments{\item{object}{A \code{pcaRes} object (result from previous PCA analysis.)}
\item{originalData}{The matrix (or ExpressionSet) that used to obtain the pcaRes object.}
\item{fold}{The number of groups to divide the data in.}
\item{nruncv}{The number of times to repeat the whole cross-validation}
\item{type}{krzanowski or imputation type cross-validation}
\item{verbose}{\code{boolean} If TRUE Q2 outputs a primitive progress bar.}
\item{Further}{arguments passed to the pca() function called within Q2.}}
\examples{data(iris)
pcIr <- pca(iris[,1:4], nPcs=3)
q2 <- Q2(pcIr, iris[,1:4], nruncv=2)
barplot(q2, xlab="Amount of PC's", ylab=expression(Q^2))}
