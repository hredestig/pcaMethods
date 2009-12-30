\name{tempFixNas}
\alias{tempFixNas}
\title{Temporary fix for missing values}
\usage{tempFixNas(mat)}
\description{Simply replace complelet missing rows or cols with zeroes.}
\value{The original matrix with completely missing rows/cols filled with zeroes.}
\author{Henning Redestig}
\arguments{\item{mat}{a matrix}}
\examples{deletediagonals(iris[,1:4], 1)}
