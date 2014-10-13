\name{leverage}
\alias{leverage}
\alias{leverage,pcaRes-method}
\title{Extract leverages of a PCA model}
\arguments{
  \item{object}{a \code{pcaRes} object}
}
\value{
  The observation leverages as a numeric vector
}
\description{
  The leverages of PCA model indicate how much influence
  each observation has on the PCA model. Observations with
  high leverage has caused the principal components to
  rotate towards them. It can be used to extract both
  "unimportant" observations as well as picking potential
  outliers.
}
\details{
  Defined as \eqn{Tr(T(T'T)^{-1}T')}{Tr(T(T'T)^(-1)T')}
}
\examples{
data(iris)
pcIr <- pca(iris[,1:4])
## versicolor has the lowest leverage
with(iris, plot(leverage(pcIr)~Species))
}
\author{
  Henning Redestig
}
\references{
  Introduction to Multi- and Megavariate Data Analysis
  using Projection Methods (PCA and PLS), L. Eriksson, E.
  Johansson, N. Kettaneh-Wold and S. Wold, Umetrics 1999,
  p. 466
}
\keyword{multivariate}

