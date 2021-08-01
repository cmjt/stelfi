#' Function to fit a self-exciting Hawkes process
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process: "mu"--base rate of the hawkes process,
#' "alpha"--intensity jump after an event occurence, and "beta"--exponential intensity decay
#' @param ... arguments to pass into \code{nlminb()}
#' @export
#' @importFrom stats optim
#' @importFrom TMB MakeADFun sdreport
setGeneric("fit_hawkes",
           function(times,parameters,...){
               standardGeneric("fit_hawkes")
           })

setMethod("fit_hawkes",
          c(times = "numeric",parameters = "vector"),
          function(times, parameters,...){
              if(!"hawkes"%in%getLoadedDLLs()){
                  dll.stelfi()
              }
              obj = MakeADFun(data = list(times = times),
                              parameters = list(log_mu = log(parameters["mu"]),
                                                logit_abratio = qlogis(parameters["alpha"] / parameters["beta"]),
                                                log_beta = log(parameters["beta"])),
                              DLL = "hawkes")
              obj$hessian <- TRUE
              opt <- optim(obj$par, obj$fn, obj$gr,...)
              return(sdreport(obj))
          })
