#' Self-exciting Hawkes process
#'
#' Fit a Hawkes process using Template Model Builder (TMB). The function \code{fit_hawkes()} fits a
#' self-exciting Hawkes process with a constant background rate. Whereas, \code{fit_hawkes_cbf()} fits a Hawkes
#' processes with a user defined custom background function (non-homogeneous background rate).
#' 
#' @details A univariate Hawkes (Hawkes, AG. 1971) process is a self-exciting temporal point process
#' with conditional intensity function
#' \eqn{\lambda(t) = \mu + \alpha \Sigma_{i:\tau_i<t}e^{(-\beta * (t-\tau_i))}}. Here \eqn{\mu} is the
#' background rate of the process and the summation describes the historic temporal dependence.  Including mark information
#' results in the conditional intensity \eqn{\lambda(t; m(t)) = \mu + \alpha \Sigma_{i:\tau_i<t}m(\tau_i)e^{(-\beta * (t-\tau_i))}},
#' where \eqn{m(t)} is the temporal mark. This model can be fitted with \code{fit_hawkes()}.
#'
#' @references Hawkes, AG. (1971) Spectra of some self-exciting and mutually exciting point processes.
#' \emph{Biometrika}, \strong{58}: 83--90.
#' 
#' @param times A vector of numeric observed time points.
#' @param parameters A named list of parameter starting values:
#' \itemize{
#' \item \code{mu}, base rate of the Hawkes process (\eqn{\mu}),
#' \item \code{alpha} (supplied if \code{model = 1}), intensity jump after an event occurrence (\eqn{\alpha}),
#' \item \code{a_par} (supplied if \code{model} = 2), logit scale for \code{alpha} (must be 
#' set so that the intensity never goes negative and so that \code{alpha} <= \code{beta}), and
#' \item \code{beta}, exponential intensity decay (\eqn{\beta}).
#' }
#' @param model A numeric indicator specifying which model to fit:
#' \itemize{
#' \item \code{model = 1}, fits a Hawkes process with exponential decay (default);
#' \item \code{model = 2}, fits a Hawkes process with an \code{alpha} that can be negative.
#' }
#' @param marks Optional, a vector of numeric marks, defaults to 1 (i.e., no marks).
#' @param tmb_silent Logical, if \code{TRUE} (default) then
#' TMB inner optimisation tracing information will be printed.
#' @param optim_silent Logical, if \code{TRUE} (default) then for each iteration
#' \code{optim()} output will be printed.
#' @param ... Additional arguments to pass to \code{optim()}
#' @return A list containing components of the fitted model, see \code{TMB::MakeADFun}. Includes
#' \itemize{
#' \item \code{par}, a numeric vector of estimated parameter values;
#' \item \code{objective}, the objective function;
#' \item \code{gr}, the TMB calculated gradient function; and
#' \item \code{simulate}, a simulation function. 
#' }
#' @examples
#' \donttest{
#' ### ********************** ###
#' ## A typical Hawkes model
#' ### ********************** ###
#' data(retweets_niwa, package = "stelfi")
#' times <- unique(sort(as.numeric(difftime(retweets_niwa, min(retweets_niwa), units = "mins"))))
#' params <- c(mu = 0.05, alpha = 0.05, beta = 0.1)
#' fit <- fit_hawkes(times = times, parameters = params)
#' get_coefs(fit)
#' ### ********************** ###
#' ## A Hawkes model with marks (ETAS-type)
#' ### ********************** ###
#' data("nz_earthquakes", package = "stelfi")
#' earthquakes <- nz_earthquakes[order(nz_earthquakes$origintime),]
#' earthquakes <- earthquakes[!duplicated(earthquakes$origintime), ]
#' times <- earthquakes$origintime
#' times <- as.numeric(difftime(times, min(times), units = "hours"))
#' marks <- earthquakes$magnitude
#' params <- c(mu = 0.05, alpha = 0.05, beta = 1)
#' fit <- fit_hawkes(times = times, parameters = params, marks = marks)
#' get_coefs(fit)
#' }
#' @seealso \code{\link{show_hawkes}}
#' @export
fit_hawkes <-  function(times, parameters = list(), model = 1,
                        marks = c(rep(1, length(times))),
                        tmb_silent = TRUE,
                        optim_silent = TRUE, ...) {
    ## general error checks
    for (i in 2:length(times)) {
        if ((times[i] - times[i - 1]) < 1.e-10)
            stop("times must be in ascending order with no simultaneous events")
    }
    if (length(marks) != length(times))
        stop("marks must have same length as times")
    if (min(marks) < 0)
        stop("marks cannot be negative")
    ## parameters
    alpha <- parameters[["alpha"]]
    if (is.null(alpha)) {
      alpha <- 0.5 * max(times) / length(times)
    }
    beta <- parameters[["beta"]]
    if (is.null(beta)) {
      beta <- 2 * alpha
    }
    mu <- parameters[["mu"]]
    if (is.null(mu)) {
      mu <- 0.5 / beta
    }
    if(model == 1) {
        ## error checks
        if (alpha > (beta / mean(marks)))
            stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
        if (alpha < 0)
            stop("alpha must be non-negative")
        ## setup
        obj <- TMB::MakeADFun(data = list(times = times, marks = marks, model_type = "hawkes"),
                              parameters = list(log_mu = log(mu),
                                                logit_abratio = stats::qlogis(alpha / (beta / mean(marks))),
                                                log_beta = log(beta)),
                              hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
    }else{
        if (model == 2) {
            ## setup
            obj <- TMB::MakeADFun(data = list(times = times, marks = marks, model_type = "neg_alpha_hawkes"),
                                  parameters = list(log_mu = log(mu),
                                                    a_par = alpha,
                                                    log_beta = log(beta)),
                                  hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
        }
    }
    trace <- if(optim_silent) 0 else 1
    opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$value
    return(obj)
}

#' @details An in-homogenous marked Hawkes process has conditional intensity function
#' \eqn{\lambda(t) = \mu(t) + \alpha \Sigma_{i:\tau_i<t}e^{(-\beta * (t-\tau_i))}}. Here, the
#' background rate, \eqn{\mu(t)}, varies in time. Such a model can be fitted
#' using \code{fit_hawkes_cbf()} where the parameters of the custom background function are estimated
#' before being passed to \code{TMB}.
#' 
#' @param background A function taking one parameter and an
#' independent variable, returning a scalar.
#' @param background_integral The integral of \code{background}.
#' @param background_parameters The parameter(s)
#' for the background function \code{background}.
#' This could be a list of multiple values.
#' @param background_min A function taking one parameter and two points,
#' returns min of \code{background} between those points.
#' @examples
#' \donttest{
#' ### ********************** ###
#' ## A Hawkes process with a custom background function
#' ### ********************** ###
#' if(require("hawkesbow")) {
#' times <- hawkesbow::hawkes(1000, fun = function(y) {1 + 0.5*sin(y)},
#' M = 1.5, repr = 0.5, family = "exp", rate = 2)$p
#' ## The background function must take a single parameter and
#' ## the time(s) at which it is evaluated
#' background <- function(params,times) {
#' A = exp(params[[1]])
#' B = stats::plogis(params[[2]]) * A
#' return(A + B  *sin(times))
#' }
#' ## The background_integral function must take a
#' ## single parameter and the time at which it is evaluated
#' background_integral <- function(params,x) {
#'         A = exp(params[[1]])
#'         B = stats::plogis(params[[2]]) * A
#'         return((A*x)-B*cos(x))
#' }
#' param = list(alpha = 0.5, beta = 1.5)
#' background_param = list(1,1)
#' fit <- fit_hawkes_cbf(times = times, parameters = param,
#' background = background,
#' background_integral = background_integral,
#' background_parameters = background_param)
#' get_coefs(fit)
#' }
#' }
#' @export
#' @rdname fit_hawkes
fit_hawkes_cbf <- function(times, parameters = list(),
                           model = 1,
                           marks = c(rep(1, length(times))),
                           background, background_integral, background_parameters,
                           background_min,
                           tmb_silent = TRUE, optim_silent = TRUE) {
    ## general error checks
    for (i in 2:length(times)) {
        if ((times[i] - times[i - 1]) < 1.e-10)
            stop("times must be in ascending order with no simultaneous events")
    }
    if (length(marks) != length(times))
        stop("marks must have same length as times")
    if (min(marks) < 0)
        stop("marks cannot be negative") 
    ## beta parameter
    beta <- parameters[["beta"]]
    if (is.null(beta)) {
      beta <- max(times) / length(times)
    } 
    if (model == 1) {
        ## alpha parameter
        alpha <- parameters[["alpha"]]
        if (is.null(alpha)) {
          alpha <- 0.5 * beta
        }
        ## error checks
        if (alpha > (beta / mean(marks)))
            stop("alpha must be smaller than or equal to beta divided by the mean of the marks")
        if (alpha < 0)
            stop("alpha must be non-negative")
        ## Nested function to be passed into optim
        optimize_background_one <- function(background_parameters, times, parameters,
                                       marks, background, background_integral,
                                       tmb_silent = TRUE, optim_silent = TRUE){
            lambda <- background(background_parameters, times)
            lambda_integral <- background_integral(background_parameters, tail(times, n = 1)) -
                background_integral(background_parameters, 0)
            
            alpha <- parameters[["alpha"]]
            beta <- parameters[["beta"]]
            obj <- TMB::MakeADFun(data = list(times = times, lambda = lambda,
                                              lambda_integral = lambda_integral,
                                              marks = marks, model_type = "custom_hawkes"),
                                  parameters = list(logit_abratio = stats::qlogis(alpha / (beta / mean(marks))),
                                                    log_beta = log(beta)),
                                  hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
            trace <- if(optim_silent) 0 else 1
            opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = 0))
            return(opt$value)
        }
        
        opt <- stats::optim(par = background_parameters, fn = optimize_background_one,
                            gr = NULL, times, parameters, marks, background,
                            background_integral,
                            tmb_silent, optim_silent)
        ## Need to run again to extract alpha and beta
        lambda <- background(opt$par, times)
        lambda_integral <- background_integral(opt$par, tail(times, n = 1)) -
            background_integral(opt$par, 0)
        alpha <- parameters[["alpha"]]
        beta <- parameters[["beta"]]
        obj <- TMB::MakeADFun(data = list(times = times, lambda = lambda,
                                          lambda_integral = lambda_integral,
                                          marks = marks, model_type = "custom_hawkes"),
                              parameters = list(logit_abratio = stats::qlogis(alpha / beta),
                                                log_beta = log(beta)),
                              hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
    }else{
        if (model == 2){
            ## a_par parameter
            a_par <- parameters[["a_par"]]
            if (is.null(a_par)) {
              a_par <- 0
            }
            ## Nested function to be passed into optim
            optimize_background_two <- function(background_parameters, times, parameters,
                                           marks, background, background_integral,
                                           background_min, tmb_silent = TRUE, optim_silent = TRUE) {
                lambda <- background(background_parameters, times)
                lambda_min <- numeric(length = length(times))
                for (k in 1:(length(times) - 1)) {
                    lambda_min[k] <- background_min(background_parameters, times[k], times[k + 1])
                }
                lambda_min[length(times)] <- background(background_parameters, tail(times, n = 1))
                if (min(lambda - lambda_min) < 0) stop("lambda_min incorrectly defined: lambda_min > lambda")
                lambda_integral <- background_integral(background_parameters, tail(times, n = 1)) -
                    background_integral(background_parameters, 0)
                
                a_par <- parameters[["a_par"]]
                beta <- parameters[["beta"]]
                obj <- TMB::MakeADFun(data = list(times = times, lambda = lambda,
                                                  lambda_min = lambda_min,
                                                  lambda_integral = lambda_integral,
                                                  marks = marks,
                                                  model_type = "neg_alpha_custom_hawkes"),
                                      parameters = list(a_par = a_par,
                                                        log_beta = log(beta)),
                                      hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
                trace <- if(optim_silent) 0 else 1
                opt <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
                return(opt$value)
            }
            
            opt <- stats::optim(par = background_parameters, fn =  optimize_background_two,
                                gr = NULL, times, parameters, marks, background,
                                background_integral,
                                background_min, tmb_silent, optim_silent)
            ## Need to run again to extract alpha and beta
            lambda <- background(opt$par, times)
            lambda_min <- numeric(length = length(times))
            for (k in 1:(length(times) - 1)) {
                lambda_min[k] <- background_min(opt$par, times[k], times[k + 1])
            }
            lambda_min[length(times)] <- background(opt$par, tail(times, n = 1))
            lambda_integral <- background_integral(opt$par, tail(times, n = 1)) -
                background_integral(opt$par, 0)
            a_par <- parameters[["a_par"]]
            beta <- parameters[["beta"]]
            ## setup
            obj <- TMB::MakeADFun(data = list(times = times, lambda = lambda,
                                              lambda_min = lambda_min,
                                              lambda_integral = lambda_integral,
                                              marks = marks,
                                              model_type = "neg_alpha_custom_hawkes"),
                                  parameters = list(a_par = a_par,
                                                    log_beta = log(beta)),
                                  hessian = TRUE, DLL = "stelfi", silent = tmb_silent)
            
        }
    }
    trace <- if(optim_silent) 0 else 1
    opt2 <- stats::optim(obj$par, obj$fn, obj$gr, control = list(trace = trace))
    obj$objective <- opt2$value
    obj$background_parameters <- opt$par
    return(obj)
}
