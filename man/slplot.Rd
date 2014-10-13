\name{slplot}
\alias{slplot}
\alias{slplot,pcaRes-method}
\title{Side by side scores and loadings plot}
\usage{
  slplot(object, pcs=c(1,2), scoresLoadings=c(TRUE, TRUE),
    sl="def", ll="def", hotelling=0.95, rug=TRUE,
    sub=NULL,...)
}
\arguments{
  \item{object}{a pcaRes object}

  \item{pcs}{which two pcs to plot}

  \item{scoresLoadings}{Which should be shown scores and or
  loadings}

  \item{sl}{labels to plot in the scores plot}

  \item{ll}{labels to plot in the loadings plot}

  \item{hotelling}{confidence interval for ellipse in the
  score plot}

  \item{rug}{logical, rug x axis in score plot or not}

  \item{sub}{Subtitle, defaults to annotate with amount of
  explained variance.}

  \item{...}{Further arguments to plot functions. Prefix
  arguments to \code{par()} with 's' for the scores plot
  and 'l' for the loadings plot. I.e. cex become scex for
  setting character expansion in the score plot and lcex
  for the loadings plot.}
}
\value{
  None, used for side effect.
}
\description{
  A common way of visualizing two principal components
}
\details{
  This method is meant to be used as a quick way to
  visualize results, if you want a more specific plot you
  probably want to get the scores, loadings with
  \code{scores(object)}, \code{loadings(object)} and then
  design your own plotting method.
}
\note{
  Uses layout instead of par to provide side-by-side so it
  works with Sweave (but can not be combined with
  \code{par(mfrow=..))}
}
\examples{
data(iris)
pcIr <- pca(iris[,1:4], scale="uv")
slplot(pcIr, sl=NULL, spch=5)
slplot(pcIr, sl=NULL, lcex=1.3, scol=as.integer(iris[,5]))
}
\author{
  Henning Redestig
}
\seealso{
  \code{\link{pca}}, \code{\link{biplot}}
}
\keyword{multivariate}

