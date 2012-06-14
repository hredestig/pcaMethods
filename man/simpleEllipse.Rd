\name{simpleEllipse}
\alias{simpleEllipse}
\title{Hotelling's T^2 Ellipse}
\usage{
  simpleEllipse(x, y, alfa = 0.95, len = 200)
}
\arguments{
  \item{x}{first variable}

  \item{y}{second variable}

  \item{alfa}{confidence level of the circle}

  \item{len}{Number of points in the circle}
}
\value{
  A matrix with X and Y coordinates for the circle
}
\description{
  Get a confidence ellipse for uncorrelated bivariate data
}
\details{
  As described in 'Introduction to multi and megavariate
  data analysis using PCA and PLS' by Eriksson et al. This
  produces very similar ellipse as compared to the ellipse
  function the ellipse package except that this function
  assumes that and y are uncorrelated (which they of are if
  they are scores or loadings from a PCA).
}
\author{
  Henning Redestig
}
\seealso{
  ellipse
}

