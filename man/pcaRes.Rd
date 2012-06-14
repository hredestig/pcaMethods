\docType{class}
\name{pcaRes}
\alias{pcaRes}
\alias{pcaRes-class}
\title{Class for representing a PCA result}
\description{
  This is a class representation of a PCA result
}
\details{
  \bold{Creating Objects}\cr \code{new("pcaRes",
  scores=[the scores], loadings=[the loadings],
  nPcs=[amount of PCs], R2cum=[cumulative R2], nObs=[amount
  of observations], nVar=[amount of variables], R2=[R2 for
  each individual PC], sDev=[stdev for each individual PC],
  centered=[was data centered], center=[original means],
  varLimit=[what variance limit was exceeded],
  method=[method used to calculate PCA], missing=[amount of
  NAs], completeObs=[estimated complete observations])}

  \bold{Slots}\cr \describe{ \item{scores}{"matrix", the
  calculated scores} \item{loadings}{"matrix", the
  calculated loadings} \item{R2cum}{"numeric", the
  cumulative R2 values} \item{sDev}{"numeric", the
  individual standard deviations of the score vectors}
  \item{R2}{"numeric", the individual R2 values}
  \item{cvstat}{"numeric", cross-validation statistics}
  \item{nObs}{"numeric", number of observations}
  \item{nVar}{"numeric", number of variables}
  \item{centered}{"logical", data was centered or not}
  \item{center}{"numeric", the original variable centers}
  \item{scaled}{"logical", data was scaled or not}
  \item{scl}{"numeric", the original variable scales}
  \item{varLimit}{"numeric", the exceeded variance limit}
  \item{nPcs,nP}{"numeric", the number of calculated PCs}
  \item{method}{"character", the method used to perform
  PCA} \item{missing}{"numeric", the total amount of
  missing values in original data}
  \item{completeObs}{"matrix", the estimated complete
  observations} \item{network}{"nlpcaNet", the network used
  by non-linear PCA} }

  \bold{Methods (not necessarily exhaustive)}\cr \describe{
  \item{print}{Print function} \item{summary}{Extract
  information about PC relevance} \item{screeplot}{Plot a
  barplot of standard deviations for PCs}
  \item{slplot}{Make a side by side score and loadings
  plot} \item{nPcs}{Get the number of PCs} \item{nObs}{Get
  the number of observations}
  \item{cvstat}{Cross-validation statistics}
  \item{nVar}{Get the number of variables}
  \item{loadings}{Get the loadings} \item{scores}{Get the
  scores} \item{dim}{Get the dimensions (number of
  observations, number of features)} \item{centered}{Get a
  logical indicating if centering was done as part of the
  model} \item{center}{Get the averages of the original
  variables.} \item{completeObs}{Get the imputed data set}
  \item{method}{Get a string naming the used PCA method}
  \item{sDev}{Get the standard deviations of the PCs}
  \item{scaled}{Get a logical indicating if scaling was
  done as part of the model} \item{scl}{Get the scales of
  the original variablesb} \item{R2cum}{Get the cumulative
  R2} }
}
\author{
  Henning Redestig
}
\keyword{classes}

