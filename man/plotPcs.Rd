\name{plotPcs}
\alias{plotPcs}
\title{Plot many side by side scores XOR loadings plots}
\usage{
  plotPcs(object, pcs = 1:nP(object),
    type = c("scores", "loadings"), sl = NULL,
    hotelling = 0.95, ...)
}
\arguments{
  \item{object}{\code{pcaRes} a pcaRes object}

  \item{pcs}{\code{numeric} which pcs to plot}

  \item{type}{\code{character} Either "scores" or
  "loadings" for scores or loadings plot respectively}

  \item{sl}{\code{character} Text labels to plot instead of
  a point, if NULL points are plotted instead of text}

  \item{hotelling}{\code{numeric} Significance level for
  the confidence ellipse. NULL means that no ellipse is
  drawn.}

  \item{...}{Further arguments to \code{\link{pairs}} on
  which this function is based.}
}
\value{
  None, used for side effect.
}
\description{
  A function that can be used to visualise many PCs plotted
  against each other
}
\details{
  Uses \code{\link{pairs}} to provide side-by-side plots.
  Note that this function only plots scores or loadings but
  not both in the same plot.
}
\examples{
data(iris)
pcIr <- pca(iris[,1:4], nPcs=3,  method="svd")
plotPcs(pcIr, col=as.integer(iris[,4]) + 1)
}
\author{
  Henning Redestig
}
\seealso{
  \code{prcomp}, \code{pca}, \code{princomp}, \code{slplot}
}
\keyword{multivariate}

