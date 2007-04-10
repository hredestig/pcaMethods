\name{pcaRes}
\docType{class}
\alias{pcaRes}
\alias{pcaRes-class}
\alias{print,pcaRes-method}
\alias{summary,pcaRes-method}
\alias{slplot,pcaRes-method}
\alias{completeObs,pcaRes-method}

\title{Class for representing a PCA result}

\description{This is a class representation of a PCA result}

\section{Creating Objects}{
  \code{new("pcaRes", scores=[the scores], loadings=[the loadings],
    nPcs=[amount of PCs], R2cum=[cumulative R2], nObs=[amount of
    observations], nVar=[amount of variables], R2=[R2 for each
    individual PC], sDev=[stdev for each individual PC],
    centered=[was data centered], center=[original means],
    varLimit=[what variance limit was exceeded], method=[method used to
    calculate PCA], missing=[amount of NAs], 
    completeObs=[estimated complete observations])}
}

\section{Slots}{
  \describe{
    \item{scores}{"matrix",  the calculated scores}
    \item{loadings}{"matrix",  the calculated loadings}
    \item{R2cum}{"numeric",  the cumulative R2 values}
    \item{sDev}{"numeric",  the individual standard
      deviations}
    \item{R2}{"numeric",  the individual R2 values}
    \item{nObs}{"numeric", amount of observations}
    \item{nVar}{"numeric", amount of variables}
    \item{centered}{"logical", data was centered or not}
    \item{center}{"numeric", the original variable centers}
    \item{varLimit}{"numeric", the exceeded variance limit}
    \item{nPcs}{"numeric", the amount of calculated PCs}
    \item{method}{"character", the method used to perform PCA}
    \item{missing}{"numeric", the total amount of missing values in
      original data}
    \item{completeObs}{"matrix", the estimated complete observations}
  }
}
\section{Methods}{
  \describe{
    \item{print}{Print function}
    \item{summary}{Extract information about PC relevance}
    \item{screeplot}{Plot a barplot of standard deviations for PCs}
    \item{slplot}{Make a side by side score and loadings plot}
    \item{biplot}{Make a scores / loadings biplot}
  }
}
\keyword{classes}
