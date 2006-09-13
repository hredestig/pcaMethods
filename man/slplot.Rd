\name{slplot}
\alias{slplot}

\title{Plot a side by side scores and loadings plot}

\description{A common way of representing PCA result for two component}
\usage{slplot(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
sl="def", ll="def", hotelling=0.95, rug=TRUE,...)}
\arguments{
  \item{object}{a pcaRes object}
  \item{pcs}{which two pcs to plot}
  \item{scoresLoadings}{Which should be shown scores and or loadings}
  \item{sl}{labels to plot in the scores  plot}
  \item{ll}{labels to plot in the loadings plot}
  \item{hotelling}{confidence interval for ellipse}
  \item{rug}{logical, rug x axis or not}
  \item{...}{Further arguments to plot functions}
}

\details{Uses layout instead of par to provide side-by-side so it works
  with Sweave.}

\value{
  None, used for side effect.
}
\author{Henning Redestig}
\seealso{\code{prcomp}, \code{pca}, \code{princomp}}
\examples{
data(iris)
pcIr <- pca(iris[,1:4], scale="UV", method="svd")
slplot(pcIr, sl=NULL, pch=5, col=as.integer(iris[,5]))
}
\keyword{multivariate}
