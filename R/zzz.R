#' Load INLA
#'
#' Loads the INLA package with \code{requireNamespace("INLA", quietly = TRUE)}.
#' On loading \code{\link{stelfi}} \code{repos['INLA']} set as
#' \code{"https://inla.r-inla-download.org/R/testing/"}.
#'
#' @references Lindgren, F. and Rue, H. (2015) Bayesian spatial modelling with R-INLA.
#' \emph{Journal of Statistical Software}, \strong{63}: 1--25.
#' @return No return value. A message is printed in the console advising if \code{INLA}
#' has been successfully loaded.
#' @export
stelfi_load_inla <- function(){
    if(requireNamespace("INLA", quietly = TRUE)) {
        message("INLA successfully loaded")
    }else{
        stop("INLA not loaded; please install using
<install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/testing'), dep=TRUE)>")
    }
}
.onLoad <- function(libname, pkgname) {
    old = options() 
    on.exit(options(old))
    repos = getOption("repos")
    repos["INLA"] = "https://inla.r-inla-download.org/R/testing"
    options(repos = repos)
    invisible(repos)
}

.onUnload <- function(libpath) {
    library.dynam.unload("stelfi", libpath)
}
