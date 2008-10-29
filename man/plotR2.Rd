\name{plotR2}
\alias{plotR2}
\title{R2 plot (screeplot) for PCA}
\description{
  Plot the R2 of the principal components to get an idea of their
  importance. Note though that the standard screeplot shows the standard
  deviations for the PC's this method shows the R2 values which empirically shows the importance of
  the PC's and is thus applicable for any PCA method rather than just
  SVD based PCA.}
\usage{
plotR2(object, nPcs=object@nPcs, type = c("barplot", "lines"), main = deparse(substitute(object)), ...)
}
\arguments{
  \item{object}{\code{pcaRes} The pcaRes object.}
  \item{nPcs}{\code{numeric} The amount of PC's to consider.}
  \item{type}{\code{character} Barplot or line plot}
  \item{main}{\code{character} The main label of the plot}
  \item{...}{Passed on to \code{screeplot}}
}
\value{
  None, used for side effect. 
}
\seealso{screeplot}
\author{Henning Redestig <redestig[at]mpimp-golm.mpg.de}
\keyword{multivariate}
