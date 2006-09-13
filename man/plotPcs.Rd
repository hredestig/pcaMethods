\name{PlotPcs}
\alias{plotPcs}

\title{Plot many side by side scores XOR loadings plots}

\description{A function that can be used to visualise many PCs plotted
  against each other}
\usage{plotPcs(object, pc=1:object@nPcs, scoresLoadings=c(TRUE, FALSE),...)}
\arguments{
  \item{object}{a pcaRes object}
  \item{pc}{which pcs to plot}
  \item{scoresLoadings}{If scores XOR loadings should be plotted}
  \item{...}{Further arguments to \code{slplot}}
}

\details{Uses par to provide side-by-side plots so it does not work
  with Sweave. Note that this function only plots scores or loadings but
  not both in the same plot.}

\value{
  None, used for side effect.
}
\author{Henning Redestig}
\seealso{\code{prcomp}, \code{pca}, \code{princomp}, \code{slplot}}
\examples{
data(iris)
pcIr <- pca(iris[,1:4], nPcs=3, scale="UV", method="svd")
plotPcs(pcIr)
}
\keyword{multivariate}
