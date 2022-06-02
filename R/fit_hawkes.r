#' Function to fit a self-exciting Hawkes process using TMB
#' @param times a vector of numeric observation time points
#' @param parameters a vector of named parmeters for the chosen model (see \code{model})
#' "mu"--base rate of the hawkes process,
#' "alpha"--intensity jump after an event occurence, and
#' "beta"--exponential intensity decay
#' @param model a factor indicator specifying which model to fit:
#' \code{1}, a hawkes process with exponential decay (default);
#' \code{2}, hawkes process with negative self-exciting intensity "jump".
#' @param marks a vector of numerical marks, defaults to 1 (i.e., no marks)
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param optim_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{optim()}
#' @export
fit_hawkes <-  function(times, parameters, model = 1,
                        marks = c(rep(1,length(times))),
                        tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    ## parameters
    mu <- parameters[["mu"]]
    alpha <- parameters[["alpha"]]
    beta <- parameters[["beta"]]
    ## general error checks
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10)
            stop("times must be in ascending order with no simultaneous events")
    }
    if (length(marks) != length(times))
        stop("marks must have same length as times")
    if (min(marks) < 0)
        stop("marks cannot be negative")
   
    if(model == 1) {
        ## error checks
        if (alpha > (beta/mean(marks)))
            stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
        if (alpha < 0)
            stop("alpha must be non-negative")
        ## check for DLL
        if (!"hawkes" %in% getLoadedDLLs()) {
            dll_stelfi("hawkes")
        }
        ## setup
        obj <- TMB::MakeADFun(data = list(times = times, marks=marks),
                              parameters = list(log_mu = log(mu),
                                                logit_abratio = stats::qlogis(alpha / (beta/mean(marks))),
                                                log_beta = log(beta)),
                              DLL = "hawkes", silent = tmb_silent)
        
    }else{
        if (model == 2){
            ## check for DLL
            if (!"neg_alpha_hawkes" %in% getLoadedDLLs()) {
                dll_stelfi("neg_alpha_hawkes")
            }
            ## setup
            obj <- TMB::MakeADFun(data = list(times = times, marks=marks),
                                     parameters = list(log_mu = log(mu),
                                                    a_par = alpha,
                                                    log_beta = log(beta)),
                                     DLL = "neg_alpha_hawkes", silent = tmb_silent)
            }
        }
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$value
    return(obj)
}

#' Function to fit a self-exciting Hawkes process with custom background function (cbf) using TMB
#' The \code{alpha} and \code{beta} parameters are estimated using TMB,
#' parameters of the CBF are optimized in R
#' @inheritParams fit_hawkes
#' @param model a factor indicator specifying which model to fit:
#' \code{1}, a hawkes process with exponential decay and cbf (default);
#' \code{2}, a hawkes process and cbf with negative self-exciting intensity "jump".
#' @param background a function taking one parameter and an independent variable, returning a scalar
#' @param background_integral the integral of background
#' @param background_parameter the parameter (which could be a list of multiple values)
#' for the background function
#' @param background_min a function taking one parameter and two points,
#' returns min of background between those points
#' @export
fit_hawkes_cbf <- function(times, parameters,
                           model = 1,
                           marks = c(rep(1,length(times))),
                           background, background_integral, background_parameters,
                           background_min, 
                           tmb_silent = TRUE, optim_silent = TRUE){
   
    ## general error checks
    for (i in 2:length(times)){
        if ((times[i]-times[i-1])<1.e-10)
            stop("times must be in ascending order with no simultaneous events")
    }
    if (length(marks) != length(times))
        stop("marks must have same length as times")
    if (min(marks) < 0)
        stop("marks cannot be negative")
    if (model == 1){
        ## parameters
        alpha <- parameters[["alpha"]]
        beta <- parameters[["beta"]]
        ## error checks
        if (alpha > (beta/mean(marks)))
            stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
        if (alpha < 0)
            stop("alpha must be non-negative")
        ## check for DLL
        if (!"custom_hawkes" %in% getLoadedDLLs()) {
            dll_stelfi("custom_hawkes")
        }
        ## Nested function to be passed into optim
        OptimizeBackground <- function(background_parameters, times, parameters,
                                       marks, background, background_integral,
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
                                  DLL = "custom_hawkes", silent = tmb_silent)
            obj$hessian <- TRUE
            trace <- if(optim_silent) 0 else 1
            opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = 0))
            return(opt$value)
        }
        
        opt <- stats::optim(par = background_parameters, fn = OptimizeBackground, 
                            gr = NULL, times, parameters, marks, background, background_integral,
                            tmb_silent, optim_silent)
        ## Need to run again to extract alpha and beta
        lambda <- background(opt$par, times)
        lambda_integral <- background_integral(opt$par, tail(times,n=1)) -
            background_integral(opt$par, 0)
        alpha <- parameters[["alpha"]]
        beta <- parameters[["beta"]]
        obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, 
                                          lambda_integral = lambda_integral, marks = marks),
                              parameters = list(logit_abratio = stats::qlogis(alpha/beta),
                                                log_beta = log(beta)),
                              DLL = "custom_hawkes", silent = tmb_silent)
    }else{
        if (model == 2){
            ## parameters
            a_par <- parameters[["a_par"]]
            beta <- parameters[["beta"]]
            ## check for DLL
            if (!"neg_alpha_custom_hawkes" %in% getLoadedDLLs()) {
                dll_stelfi("neg_alpha_custom_hawkes")
            }
                                        # Nested function to be passed into optim
            OptimizeBackground <- function(background_parameters, times, parameters,
                                           marks, background, background_integral,
                                           background_min, tmb_silent = TRUE, optim_silent = TRUE){
                lambda <- background(background_parameters, times)
                lambda_min <- numeric(length=length(times))
                for (k in 1:(length(times)-1)){
                    lambda_min[k] <- background_min(background_parameters,times[k],times[k+1])
                }
                lambda_min[length(times)] <- background(background_parameters, tail(times, n=1))
                if (min(lambda-lambda_min) < 0) stop("lambda_min incorrectly defined: lambda_min > lambda")
                lambda_integral <- background_integral(background_parameters, tail(times,n=1)) -
                    background_integral(background_parameters, 0)
                
                a_par <- parameters[["a_par"]]
                beta <- parameters[["beta"]]
                obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, lambda_min = lambda_min,
                                                  lambda_integral = lambda_integral, marks = marks),
                                      parameters = list(a_par = a_par,
                                                        log_beta = log(beta)),
                                      DLL = "neg_alpha_custom_hawkes", silent = tmb_silent)
                obj$hessian <- TRUE
                trace <- if(optim_silent) 0 else 1
                opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
                return(opt$value)
            }
            
            opt <- stats::optim(par = background_parameters, fn = OptimizeBackground, 
                                gr = NULL, times, parameters, marks, background, background_integral,
                                background_min, tmb_silent, optim_silent)
                                        # Need to run again to extract alpha and beta
            lambda <- background(opt$par, times)
            lambda_min <- numeric(length=length(times))
            for (k in 1:(length(times)-1)){
                lambda_min[k] <- background_min(opt$par,times[k],times[k+1])
            }
            lambda_min[length(times)] <- background(opt$par, tail(times, n=1))
            lambda_integral <- background_integral(opt$par, tail(times,n=1)) -
                background_integral(opt$par, 0)
            a_par <- parameters[["a_par"]]
            beta <- parameters[["beta"]]
            ## setup
            obj <- TMB::MakeADFun(data = list(times = times,lambda = lambda, lambda_min = lambda_min,
                                              lambda_integral = lambda_integral, marks = marks),
                                  parameters = list(a_par = a_par,
                                                    log_beta = log(beta)),
                                  DLL = "neg_alpha_custom_hawkes", silent = tmb_silent)
            
        }
    }
    obj$hessian <- TRUE
    trace <- if(optim_silent) 0 else 1
    opt2 <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
    res <- TMB::sdreport(obj)
    return(list(alpha = res$value[1], beta = res$value[2], background_parameters = opt$par,
                objective = opt2$value))
}


