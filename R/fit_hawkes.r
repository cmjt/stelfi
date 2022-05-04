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
    # Verify arguments
    mu <- parameters[["mu"]]
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    
    if (alpha > beta) stop("alpha must be smaller than or equal to beta")
    if (alpha < 0) stop("alpha must be non-negative")
    times_sorted <- sort(times)
    if (!setequal(times,times_sorted)) stop("times must be in ascending order")
    
    time_elapsed <- numeric(length=length(times)-1)
    for (i in 1:length(time_elapsed)){
        time_elapsed[i] = times[i+1]-times[i]
    }
    if (min(time_elapsed) == 0) stop("Concurrent events not permitted in Poisson processes")
    
    if (!"hawkes" %in% getLoadedDLLs()) {
        dll_stelfi("hawkes")
    }
    obj <- TMB::MakeADFun(data = list(times = times),
                          parameters = list(log_mu = log(mu),
                                            logit_abratio = stats::qlogis(alpha / beta),
                                            log_beta = log(beta)),
                          DLL = "hawkes", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    return(obj)
}
