\name{R2VX}
\alias{R2VX}
\alias{R2VX,pcaRes-method}
\title{R2 goodness of fit}
\arguments{
  \item{object}{a PCA model object}

  \item{direction}{choose between calculating R2 per
  variable, per observation or for the entire data with
  'variables', 'observations' or 'complete'.}

  \item{data}{the data used to fit the model}

  \item{pcs}{the number of PCs to use to calculate R2}
}
\value{
  A vector with R2 values
}
\description{
  Flexible calculation of R2 goodness of fit.
}
\examples{
R2VX(pca(iris))
}
\author{
  Henning Redestig
}

