\name{plotR2}
\alias{plotR2}
\title{R2 plot (screeplot) for PCA}
\usage{plotR2(object, nPcs=nP(object), type=c("barplot", "lines"),
    main=deparse(substitute(object)), ...)}
\description{Plot the R2 of the principal components to get an idea of their
importance. Note though that the standard screeplot shows the
standard deviations for the PC's this method shows the R2 values
which empirically shows the importance of the PC's and is thus
applicable for any PCA method rather than just SVD based PCA.}
\note{This method is deprecated in favor of plot.pcaRes which does
(almost) the same thing but with a better name.}
\value{None, used for side effect.}
\seealso{\link{screeplot}}
\keyword{multivariate}
\author{Henning Redestig}
\arguments{\item{object}{\code{pcaRes} The pcaRes object.}
\item{nPcs}{\code{numeric} The number of PC's to consider.}
\item{type}{\code{character} Barplot or line plot}
\item{main}{\code{character} The main label of the plot}
\item{...}{Passed on to \code{screeplot}}}
