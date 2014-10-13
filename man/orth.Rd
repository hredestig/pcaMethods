\name{orth}
\alias{orth}
\title{Calculate an orthonormal basis}
\usage{
  orth(mat, skipInac = FALSE)
}
\arguments{
  \item{mat}{matrix to calculate orthonormal base}

  \item{skipInac}{do not include components with precision
  below .Machine$double.eps if TRUE}
}
\value{
  orthonormal basis for the range of matrix
}
\description{
  ONB = orth(mat) is an orthonormal basis for the range of
  matrix mat.  That is, ONB' * ONB = I, the columns of ONB
  span the same space as the columns of mat, and the number
  of columns of ONB is the rank of mat.
}
\author{
  Wolfram Stacklies
}

