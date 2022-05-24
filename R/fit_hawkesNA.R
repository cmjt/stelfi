#' Function to fit a self-exciting Hawkes process using TMB
#' Alpha can be negative
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process:
#' "mu"--base rate of the hawkes process,
#' "a_par"--controls alpha, see TMB template
#' "beta"--exponential intensity decay
#' @param marks a vector of numerical marks, defaults to ones
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{optim()}
#' @export
fit_hawkesNA <-  function(times, parameters, marks=c(rep(1,length(times))), tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    # Verify arguments
    mu <- parameters[["mu"]]
    a_par <- parameters[["a_par"]]
    beta <- parameters[["beta"]]
    
    
    
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10) stop("times must be in ascending order with no simultaneous events")
    }
    
    if (length(marks) != length(times)) stop("marks must have same length as times")
    
    if (!"hawkesNA" %in% getLoadedDLLs()) {
        dll_stelfi("hawkesNA")
    }
    obj <- TMB::MakeADFun(data = list(times = times, marks=marks),
                          parameters = list(log_mu = log(mu),
                                            a_par = a_par,
                                            log_beta = log(beta)),
                          DLL = "hawkesNA", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$value
    return(obj)
}

#' Function to fit a self-exciting Hawkes process with custom background function (CBF) using TMB
#' The alpha and beta parameters are minimized with TMB, parameters of the CBF are optimized in R
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process:
#' "a_par"--controls alpha, see TMB template
#' "beta"--exponential intensity decay
#' @param marks a vector of numerical marks, defaults to ones
#' @param background a function taking one parameter and an independent variable, returning a scalar
#' @param background_integral the integral of background
#' @param background_min a function taking one parameter and two points, returns min of background between those points
#' @param background_parameter the parameter (which could be a list of multiple values) for the background function
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @export
#' 
fit_hawkesNACBF <- function(times, parameters, marks=c(rep(1,length(times))),
                          background, background_integral, background_min, background_parameters,
                          tmb_silent = TRUE, optim_silent = TRUE){
    # Argument verification
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10) stop("times must be in ascending order with no simultaneous events")
    }
    a_par <- parameters[["a_par"]]
    beta <- parameters[["beta"]]
    #if (alpha > beta) stop("alpha must be smaller than or equal to beta")
    #if (alpha < 0) stop("alpha must be non-negative")
    
    if (length(marks) != length(times)) stop("marks must have same length as times")
    
    if (!"hawkesNACBF" %in% getLoadedDLLs()) {
        dll_stelfi("hawkesNACBF")
    }
    # Nested function to be passed into optim
    OptimizeBackground <- function(background_parameters, times, parameters, marks, background, background_integral,
                                   background_min, tmb_silent = TRUE, optim_silent = TRUE){
        lambda <- background(background_parameters, times)
        lambda_min <- numeric(length=length(times))
        for (k in 1:(length(times)-1)){
            lambda_min[k] <- background_min(background_parameters,times[k],times[k+1])
        }
        lambda_min[length(times)] <- background_min(background_parameters, times[length(times)],times[length(times)])
        if (min(lambda-lambda_min) < 0) stop("lambda_min incorrectly defined: lambda_min[i] > lambda[i]")
        lambda_integral <- background_integral(background_parameters, tail(times,n=1)) -
            background_integral(background_parameters, 0)
        
        a_par <- parameters[["a_par"]]
        beta <- parameters[["beta"]]
        obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, lambda_min = lambda_min,
                                          lambda_integral = lambda_integral, marks = marks),
                              parameters = list(a_par = a_par,
                                                log_beta = log(beta)),
                              DLL = "hawkesNACBF", silent = tmb_silent)
        obj$hessian <- TRUE
        trace <- if(optim_silent) 0 else 1
        opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
        return(opt$value)
    }
    
    opt <- stats::optim(par = background_parameters, fn = OptimizeBackground, 
                        gr = NULL, times, parameters, marks, background, background_integral,
                        background_min, tmb_silent, optim_silent)
    # Need to run again to extract alpha and beta
    print(opt$par)
    print(opt$value)
    lambda <- background(opt$par, times)
    lambda_min <- numeric(length=length(times))
    for (k in 1:(length(times)-1)){
        lambda_min[k] <- background_min(opt$par,times[k],times[k+1])
    }
    lambda_min[length(times)] <- background_min(opt$par, times[length(times)],times[length(times)])
    lambda_integral <- background_integral(opt$par, tail(times,n=1)) -
        background_integral(opt$par, 0)
    a_par <- parameters[["a_par"]]
    beta <- parameters[["beta"]]
    obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, lambda_min = lambda_min,
                                      lambda_integral = lambda_integral, marks = marks),
                          parameters = list(a_par = a_par,
                                            log_beta = log(beta)),
                          DLL = "hawkesNACBF", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt2 <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
    res <- TMB::sdreport(obj)
    return(list(alpha=res$value[1],beta=res$value[2],background_parameters=opt$par,
                obj=opt2$value))
}