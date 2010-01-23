\name{DModX}
\title{DModX}
\usage{DModX(object, dat, ...)}
\description{Distance to the model of X-space.}
\details{Measures how well described the observations are, i.e. how well
they fit in the mode. High DModX indicate a poor fit. Defined as:

\eqn{\frac{\sqrt{\frac{SSE_i}{K-A}}}{\sqrt{\frac{SSE}{(N-A-A_0)(K-A)}}}}

For observation \eqn{i}, in a model with \eqn{A} components,
\eqn{K} variables and \eqn{N} obserations. SSE is the squared sum
of the residuals. \eqn{A_0} is 0 if model was centered and 1
otherwise. DModX is claimed to be approximately F-distributed and
can therefore be used to check if an observation is significantly
far away from the PCA model assuming normally distributed data.

Pass original data as an argument if the model was calculated with
\code{completeObs=FALSE}.}
\value{A vector with distances from observations to the PCA model}
\alias{DModX}
\alias{DModX,pcaRes-method}
\references{Introduction to Multi- and Megavariate Data Analysis
using Projection Methods (PCA and PLS), L. Eriksson, E. Johansson,
N. Kettaneh-Wold and S. Wold, Umetrics 1999, p. 468}
\author{Henning Redestig}
\arguments{\item{object}{a pcaRes object}
\item{dat}{the original data, taken from \code{completeObs} if
left missing.}
\item{...}{Not used}}
\examples{data(iris)
pcIr <- pca(iris[,1:4])
with(iris, plot(DModX(pcIr)~Species))}
