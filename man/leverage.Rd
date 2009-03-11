\name{leverage}
\alias{leverage}

\title{Extract leverages of a PCA model}

\description{The leverages of PCA model indicate how much influence each
observation has on the PCA model. Observations with high leverage has
caused the principal components to rotate towards them. It can be used
to extract both "unimportant" observations as well as picking potential
outliers.}

\usage{leverage(object,...)}
\arguments{
  \item{object}{a pcaRes object}
  \item{...}{not used}
}

\details{Defined as \eqn{Tr(T(T'T)^{-1}T')}{Tr(T(T'T)^(-1)T')}}

\value{
  The observation leverages as a numeric vector
}
\references{Introduction to Mult- and Megavaraite Data Analysis uing
  Projection Methods (PCA and PLS), L. Eriksson, E. Johansson,
  N. Kettaneh-Wold and S. Wold, Umetrics 1999, p. 466}
\author{Henning Redestig}
\examples{
data(iris)
pcIr <- pca(iris[,1:4])
## versicolor has the lowest leverage
plot(leverage(pcIr)~iris$Species)
}
\keyword{multivariate}
