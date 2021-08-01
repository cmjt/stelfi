#' creates a virtual class that is a superclass to the component classes so then both children inherit from that class
setClassUnion("missingOrNULL", c("missing", "NULL")) 
#' Fuction to compile TMB C++ templates
#'
#' Compiles the  TMB templates into a shared object file.
#' Must be done a single time following installation or updating of the package.
#' @export
setGeneric("compile.stelfi",
           function(x = NULL){
               standardGeneric("compile.stelfi")
           })
setMethod("compile.stelfi",
          c(x = "missingOrNULL"),
          function(x){
              wd <- getwd()
              dir <- paste(system.file(package = "stelfi"), "/src", sep = "")
              setwd(dir)
              if (!dir.exists("../bin")){
                  dir.create("../bin")
              }
              files <- strsplit(list.files(), "[.]")
              base <- sapply(files, function(x) x[1])
              ext <- sapply(files, function(x) x[2]) 
              is.windows <- length(grep("Windows", sessionInfo()$running)) > 0
              for (i in base[ext == "cpp"]){
                  compile(paste(i, ".cpp", sep = ""))
                  unlink(paste(i, ".o", sep = ""))
                  if (is.windows)
                      file.rename(paste(i, ".dll", sep = ""),
                              paste("../bin/", i, ".dll", sep = ""))
                  else
                      file.rename(paste(i, ".so", sep = ""),
                                  paste("../bin/", i, ".so", sep = ""))
              }
              setwd(wd)
          })

#' Function to load DLLs for C++ templtes
#' Loads required DLLs for models fitted using TMB
#' @export
setGeneric("dll.stelfi",
           function(x = NULL){
               standardGeneric("dll.stelfi")
           })
setMethod("dll.stelfi",
          c(x = "missingOrNULL"),
          function(x){
              dll.dir <- paste(system.file(package = "stelfi"), "/bin/", sep = "")
              for (i in paste(dll.dir, list.files(dll.dir), sep = "")){
                  dyn.load(i)
              }
          })
