.onLoad <- function(lib, pkg) {
    cat("Loading compiled code...\n")
    library.dynam("stelfi", pkg, lib)
}
.onUnload <- function(libpath) {
  library.dynam.unload("stelfi", libpath)
}
