#' Fuction to compile TMB C++ templates
#'
#' Compiles the  TMB templates into a shared object file.
#' Must be done a single time following installation or updating of the package.
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
        if (is_windows)
            file.rename(paste(i, ".dll", sep = ""),
                        paste("../bin/", i, ".dll", sep = ""))
        else
            file.rename(paste(i, ".so", sep = ""),
                        paste("../bin/", i, ".so", sep = ""))
    }
    setwd(wd)
}

#' Function to load DLLs for C++ templtes
#' Loads required DLLs for models fitted using TMB
#' @param x option, if provided specifies \code{stelfi} DLL to load
#' @export
dll_stelfi <- function(x) {
    dll_dir <- paste(system.file(package = "stelfi"), "/bin/", sep = "")
    if(missing(x)){
        for (i in paste(dll_dir, list.files(dll_dir), sep = "")) {
            dyn.load(i)
        }
    }else{
        dyn.load(paste(dll_dir, x, sep = "/"))
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
