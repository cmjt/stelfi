#' Function to fit a self-exciting Hawkes process
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process: "mu"--base rate of the hawkes process,
#' "alpha"--intensity jump after an event occurence, and "beta"--exponential intensity decay
#' @param ... arguments to pass into nlminb
#' @export
#'
setGeneric("fit.hawkes",
           function(times,parameters,...){
               standardGeneric("fit.hawkes")
           })

setMethod("fit.hawkes",
          c(times = "numeric",parameters = "vector"),
          function(times, parameters,...){
              if(!"hawkes"%in%getLoadedDLLs()){
                  dll.stelfi()
              }
              obj = TMB::MakeADFun(data = list(times = times), parameters = parameters,DLL = "hawkes")
              opt = stats::nlminb(obj$par,obj$fn,obj$gr,...)
              return(obj)
          })
