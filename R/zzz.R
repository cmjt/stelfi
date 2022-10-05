#' Load INLA package
#'
#' Loads the INLA package with \code{requireNamespace("INLA", quietly = TRUE)}.
#' On loading \code{\link{stelfi}} \code{repos['INLA']} set as
#' \url{"https://inla.r-inla-download.org/R/testing"}.
#'
#' @export
stelfi_load_inla <- function(){
    if(requireNamespace("INLA", quietly = TRUE)) {
        message("INLA successfully loaded")
    }else{
        message("INLA not loaded; please install using
install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/testing'), dep=TRUE)")
    }
}
.onLoad <- function(libname) {
  repos = getOption("repos")
  repos["INLA"] = "https://inla.r-inla-download.org/R/testing"
  options(repos = repos)
  invisible(repos)
}

.onUnload <- function(libpath) {
  library.dynam.unload("stelfi", libpath)
}
