#' Function to compile all TMB C++ templates contained in the package
#' 
#' \code{compile_stelfi} compiles the  TMB templates into a shared object file.
#' Must be done a single time following installation or updating of the package. For
#' loading the required DLL see \code{\link{dll_stelfi}}.
#' 
#' @seealso \code{\link{dll_stelfi}}
#' @examples \dontrun{
#' compile_stelfi()
#' }
#' @export
compile_stelfi <- function() {
    wd <- getwd()
    dir <- paste(system.file(package = "stelfi"), "/src", sep = "")
    setwd(dir)
    if (!dir.exists("../bin")) {
        dir.create("../bin")
    }
    files <- strsplit(list.files(), "[.]")
    base <- sapply(files, function(x) x[1])
    ext <- sapply(files, function(x) x[2])
    is_windows <- length(grep("Windows", utils::sessionInfo()$running)) > 0
    for (i in base[ext == "cpp"]) {
        TMB::compile(paste(i, ".cpp", sep = ""))
        unlink(paste(i, ".o", sep = ""))
        if (is_windows) {
            file.rename(paste(i, ".dll", sep = ""),
                        paste("../bin/", i, ".dll", sep = ""))
        } else {
            file.rename(paste(i, ".so", sep = ""),
                        paste("../bin/", i, ".so", sep = ""))
            }
    }
    setwd(wd)
}

#' Function to load DLLs for C++ templates
#' 
#' \code{dll_stelfi} loads required DLLs for models fitted using TMB.
#' For compiling the required templates see \code{\link{compile_stelfi}}.
#' 
#' @param x Optional, if provided specifies which \code{\link{stelfi}} DLL to load.
#' @examples \dontrun{
#' dll_stelfi()
#' }
#' @export
dll_stelfi <- function(x) {
    dll_dir <- paste(system.file(package = "stelfi"), "/bin/", sep = "")
    is_windows <- length(grep("Windows", utils::sessionInfo()$running)) > 0
    if(missing(x)){
        for (i in paste(dll_dir, list.files(dll_dir), sep = "")) {
            dyn.load(i)
        }
    } else if (is_windows) {
        dyn.load(paste(dll_dir, x, ".dll", sep = ""))
    } else {
        dyn.load(paste(dll_dir, x, ".so", sep = ""))
    }
}
#' @importFrom TMB compile MakeADFun sdreport
#' @import Rcpp
#' @import Matrix
#' @importFrom ggplot2 ggplot aes xlab geom_line theme_minimal geom_histogram
#' @importFrom gridExtra grid.arrange
#' @import rgeos
#' @import sf
#' @import sp
#' @import spatstat.core
#' @import spatstat.geom
#' @importFrom stats optim qlogis runif
#' @importFrom utils sessionInfo
#' @importFrom grDevices dev.new
#' @import INLA
#' @importFrom grDevices dev.new
#' @importFrom graphics points
#' @importFrom stats ppoints qexp quantile runif
#' @importFrom utils tail
