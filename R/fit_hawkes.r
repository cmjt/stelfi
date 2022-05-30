#' Function to fit a self-exciting Hawkes process using TMB
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the hawkes process:
#' "mu"--base rate of the hawkes process,
#' "alpha"--intensity jump after an event occurence, and
#' "beta"--exponential intensity decay
#' @param marks a vector of numerical marks, defaults to ones
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{optim()}
#' @export
fit_hawkes <-  function(times, parameters, marks=c(rep(1,length(times))), tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    # Verify arguments
    mu <- parameters[["mu"]]
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    
    if (alpha > (beta/mean(marks))) stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
    if (alpha < 0) stop("alpha must be non-negative")
    
    
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10) stop("times must be in ascending order with no simultaneous events")
    }
    
    if (length(marks) != length(times)) stop("marks must have same length as times")
    if (min(marks) < 0) stop("marks cannot be negative")
    
    if (!"hawkes" %in% getLoadedDLLs()) {
        dll_stelfi("hawkes")
    }
    obj <- TMB::MakeADFun(data = list(times = times, marks=marks),
                          parameters = list(log_mu = log(mu),
                                            logit_abratio = stats::qlogis(alpha / (beta/mean(marks))),
                                            log_beta = log(beta)),
                          DLL = "hawkes", silent = tmb_silent)
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
#' "alpha"--intensity jump after an event occurence, and
#' "beta"--exponential intensity decay
#' @param marks a vector of numerical marks, defaults to ones
#' @param background a function taking one parameter and an independent variable, returning a scalar
#' @param background_integral the integral of background
#' @param background_parameter the parameter (which could be a list of multiple values) for the background function
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @export
#' 
fit_hawkesCBF <- function(times, parameters, marks=c(rep(1,length(times))),
                          background, background_integral, background_parameters,
                          tmb_silent = TRUE, optim_silent = TRUE){
    # Argument verification
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10) stop("times must be in ascending order with no simultaneous events")
    }
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    if (alpha > (beta/mean(marks))) stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
    if (alpha < 0) stop("alpha must be non-negative")
    
    if (length(marks) != length(times)) stop("marks must have same length as times")
    if (min(marks) < 0) stop("marks cannot be negative")
    
    if (!"hawkesCBF" %in% getLoadedDLLs()) {
        dll_stelfi("hawkesCBF")
    }
    # Nested function to be passed into optim
    OptimizeBackground <- function(background_parameters, times, parameters, marks, background, background_integral,
                                   tmb_silent = TRUE, optim_silent = TRUE){
        lambda <- background(background_parameters, times)
        lambda_integral <- background_integral(background_parameters, tail(times,n=1)) -
            background_integral(background_parameters, 0)
        
        alpha <- parameters[["alpha"]]
        beta <- parameters[["beta"]]
        obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, 
                                          lambda_integral = lambda_integral, marks = marks),
                              parameters = list(logit_abratio = stats::qlogis(alpha/(beta/mean(marks))),
                                                log_beta = log(beta)),
                              DLL = "hawkesCBF", silent = tmb_silent)
        obj$hessian <- TRUE
        trace <- if(optim_silent) 0 else 1
        opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = 0))
        return(opt$value)
    }
    
    opt <- stats::optim(par = background_parameters, fn = OptimizeBackground, 
                        gr = NULL, times, parameters, marks, background, background_integral,
                        tmb_silent, optim_silent)
    # Need to run again to extract alpha and beta
    lambda <- background(opt$par, times)
    lambda_integral <- background_integral(opt$par, tail(times,n=1)) -
        background_integral(opt$par, 0)
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, 
                                      lambda_integral = lambda_integral, marks = marks),
                          parameters = list(logit_abratio = stats::qlogis(alpha/beta),
                                            log_beta = log(beta)),
                          DLL = "hawkesCBF", silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt2 <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
    res <- TMB::sdreport(obj)
    return(list(alpha=res$value[1],beta=res$value[2],background_parameters=opt$par,
                objective=opt2$value))
}

