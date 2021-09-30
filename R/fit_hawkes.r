#' Function to fit a self-exciting Hawkes process using TMB
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process:
#' "mu"--base rate of the hawkes process,
#' "alpha"--intensity jump after an event occurence, and
#' "beta"--exponential intensity decay
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{optim()}
#' @export
fit_hawkes <-  function(times, parameters, tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    if (!"hawkes" %in% getLoadedDLLs()) {
        dll_stelfi("hawkes")
    }
    obj <- TMB::MakeADFun(data = list(times = times),
                          parameters = list(log_mu = log(parameters["mu"]),
                                            logit_abratio = stats::qlogis(parameters["alpha"] / parameters["beta"]),
                                            log_beta = log(parameters["beta"])),
                          DLL = "hawkes", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    return(obj)
}
