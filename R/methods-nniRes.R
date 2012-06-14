##' Print a brief description of nniRes model
##' @title Print a nniRes model
##' @param x An \code{nniRes} object
##' @param ... Not used
##' @return Nothing, used for side-effect
##' @export
##' @author Henning Redestig
showNniRes <- function(x, ...) {
  summary(x)
  cat(dim(x)["nVar"], "\tVariables\n")
  cat(dim(x)["nObs"],"\tSamples\n")
  cat(nmissing(x), "\tNAs (",
      round(100 * nmissing(x) / (nObs(x) * nVar(x)),
            getOption("str")$digits.d), "%)\n")
  cat("k was set to", x@k, "\n")
  if(centered(x))
    cat("Data was mean centered before running LLSimpute \n")
  else
    cat("Data was NOT mean centered before running LLSimpute \n")
  if(scaled(x))
    cat("Data was scaled before running LLSimpute \n")
  else
    cat("Data was NOT scaled before running LLSimpute \n")
}
setMethod("print", "nniRes", showNniRes)
setMethod("show", "nniRes", function(object) showNniRes(object))
