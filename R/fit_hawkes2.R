#' Function to fit a self-exciting Hawkes process using TMB with non-constant background
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
fit_hawkes2 <-  function(times, parameters, tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    # Verify arguments
    mu <- parameters[["mu"]]
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    a <- parameters[["a"]]
    b <- parameters[["b"]]
    c <- parameters[["c"]]
    
    if (alpha > beta) stop("alpha must be smaller than or equal to beta")
    if (alpha < 0) stop("alpha must be non-negative")
    if ((c > pi) || (c < -pi)) stop("The absolute value of C must be less than or equal to pi")
    if (a > mu) stop("a must be smaller than or equal to mu")
    if (a < 0) stop("a must be non-negative")
    
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10) stop("times must be in ascending order with no simultaneous events")
    }
    
    if (!"hawkes2" %in% getLoadedDLLs()) {
        dll_stelfi("hawkes2")
    }
    obj <- TMB::MakeADFun(data = list(times = times),
                          parameters = list(log_mu = log(mu),
                                            logit_abratio = stats::qlogis(alpha/beta),
                                            log_beta = log(beta),
                                            logit_amuratio = stats::qlogis(a/mu),
                                            b = b,
                                            atanh_c = atanh(c/pi)),
                          DLL = "hawkes2", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$value
    return(obj)
}
