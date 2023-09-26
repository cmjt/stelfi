#' Unload stelfi dlls
.onUnload <- function(libpath) {
    library.dynam.unload("stelfi", libpath)
}
