#' Function to fit a spde spatial Hawkes process using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' @param coefs the logged base rate of the Hawkes process and coefficients of covariates
#' @param designmat design matrix
#' @param logit_abratio logit ration of a, b parameters of the temporal Hawkes process
#' @param log_beta the \code{log(beta)} base rate of the Hawkes process
#' @param log_kappa the \code{log(kappa)} parameter for the GMRF
#' @param log_tau the \code{log(tau)} parameter for the GMRF
#' @param log_xsigma the log standard deviation on x-axis of self-exciting kernel
#' (if \code{time_independent = FALSE}, it is the s.d. after 1 time unit)
#' @param log_ysigma the log  standard deviation on y-axis of self-exciting kernel
#' (if \code{time_independent = FALSE}, it is the s.d. after 1 time unit)
#' @param tmax maximum time (\code{numeric})
#' @param reltol \code{numeric}, relative tolerance (default  \code{1e-12})
#' @param abstol \code{numeric}, absolute tolerance (default  \code{1e-12})
#' @param lmat \code{INLA::inla.spde.make.A(smesh, locs)}
#' @param simple  binary, should Gaussian kernels have a
#' covariate matrix that is proportional to time since the event.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
#' @inheritParams fit_lgcp_tmb
#' @noRd
fit_hspde_tmb <- function(times, locs, sf, smesh,
                          coefs, designmat, logit_abratio = 0, log_beta = 0,
                          log_kappa = 0, log_tau = 0, log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0,
                          spde, w, tmax = max(times),
                          reltol = 1e-12, abstol = 1e-12,
                          lmat = lmat, simple = simple,
                          tmb_silent, nlminb_silent, ...) {
    innerloc <- Reduce(rbind, apply(smesh$graph$tv, 1, function(vt) {
        temp <-  sf::st_polygon(list(smesh$loc[rep(vt, length = length(vt) + 1), 1:2]))
        if (is.null(sf::st_intersection(temp, sf)))
          rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs,
                xyloc = smesh$loc[, 1:2], reltol = reltol, abstol = abstol,
                spde = spde$param.inla[c("M0", "M1", "M2")], w = w, tmax = tmax,
                designmat = designmat, tv = innerloc, simple = simple, lmat = lmat,
                model_type = "spde_hawkes")
    param <- list(coefs = coefs, logit_abratio = logit_abratio, log_beta = log_beta,
                  log_xsigma =  log_xsigma, log_ysigma =  log_ysigma,
                  atanh_rho = atanh_rho, log_kappa = log_kappa, log_tau = log_tau,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = 1))
    obj <- TMB::MakeADFun(data, param, hessian = TRUE, random = c("x"),
                           DLL = "stelfi", silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$objective
    return(obj)
}
#' Function to fit a spatial Hawkes process using \code{TMB}
#' @inheritParams fit_hspde_tmb
#' @noRd
fit_hspat_tmb <- function(times, locs, sf,
                          smesh, coefs, designmat, logit_abratio = 0, log_beta = 0,
                          log_xsigma = 0, log_ysigma = 0, atanh_rho = 0, w,
                          reltol = 1e-12, abstol = 1e-12,
                          tmax = max(times),
                          lmat, simple,
                          tmb_silent, nlminb_silent, ...) {
    innerloc <- Reduce(rbind, apply(smesh$graph$tv, 1, function(vt) {
        temp <-  sf::st_polygon(list(smesh$loc[rep(vt, length = length(vt) + 1), 1:2]))
        if (is.null(sf::st_intersection(temp, sf)))
          rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs, 
                 xyloc = smesh$loc[,1:2], reltol = reltol, abstol = abstol,
                  w = w, tmax = tmax, designmat = designmat,
                 tv = innerloc, simple = simple, lmat = lmat,
                 model_type = "spatial_hawkes")
    param <- list(coefs = coefs, logit_abratio = logit_abratio,
                  log_beta = log_beta,
                  log_xsigma =  log_xsigma,
                  log_ysigma =  log_ysigma, atanh_rho = atanh_rho)
    obj <- TMB::MakeADFun(data, param, hessian = TRUE,
                           DLL = "stelfi", silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$objective
    return(obj)
}
#' Modelling spatiotemporal self-excitement
#' 
#' Fits spatiotemporal Hawkes models. The self-excitement is 
#' Gaussian in space and exponentially decaying in time.
#'
#' @param locs A \code{data.frame} of \code{x} and \code{y} locations, 2xn.
#' @param parameters A list of named parameters:
#' \itemize{
#' \item \code{coefs}, logged base rate of the Hawkes process and coefficients of covariates
#' \item \code{alpha}, intensity jump after an event occurrence
#' \item \code{beta},  rate of exponential decay of intensity after event occurrence
#' \item \code{tau},  \eqn{\tau} parameter for the GMRF (supplied only if \code{GMRF = TRUE})
#' \item \code{kappa}, \eqn{\kappa} parameter for the GMRF (supplied only if \code{GMRF = TRUE})
#' \item \code{xsigma}, standard deviation on x-axis of self-exciting kernel (if \code{time_independent = FALSE}, it is the s.d. after 1 time unit)
#' \item \code{ysigma}, standard deviation on y-axis of self-exciting kernel (if \code{time_independent = FALSE}, it is the s.d. after 1 time unit)
#' \item \code{rho}, correlation between x and y for the self-exciting kernel (the off-diagonal elements in the kernel's covariate matrix are \code{xsigma * ysigma * rho})
#' }
#' @param covariates Optional, a \code{matrix} of covariates at each
#' \code{smesh} node.
#' @param GMRF Logical, default \code{FALSE}. If \code{TRUE}, a Gaussian Markov
#'  Random Field is included as a latent spatial effect.
#' @param time_independent Logical, default \code{TRUE}. If \code{FALSE}, Gaussian kernels have a
#' covariate matrix that is proportional to time since the event.
#' Warning, this is very memory intensive.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
#' 
#' @details Temporal self-excitement follows an exponential decay function.
#' The self-excitement over space follows a Gaussian distribution centered at the triggering event.
#' There are two formulations of this model. The default is that the Gaussian function has a fixed spatial
#' covariance matrix, independent of time. Alternatively, covariance can be directly proportional to time,
#' meaning that the self-excitement radiates out from the center over time.
#' This can be appropriate when the mechanism causing self-excitement travels
#' at a finite speed, but is very memory-intensive. The spatiotemporal intensity function
#' used by \code{\link{stelfi}} is
#' \eqn{\lambda(s,t) = \mu + \alpha \Sigma_{i:\tau_i<t}(\textrm{exp}(-\beta * (t-\tau_i)) G_i(s-x_i, t - \tau_i))}
#' where
#' \itemize{
#' \item \eqn{\mu} is the background rate,
#' \item \eqn{\beta} is the rate of temporal decay,
#' \item \eqn{\alpha} is the increase in intensity after an event,
#' \item \eqn{\tau_i} are the event times,
#' \item \eqn{x_i} are the event locations (in 2D Euclidean space), and
#' \item \eqn{G_i(s-x_i, t - \tau_i)} is the spatial self-excitement kernel.
#' }
#' \eqn{G_i(.,.)} can take two forms:
#' \itemize{
#' \item For time-independent spatial excitement (\code{time_independent = TRUE}),
#' \eqn{G_i(s-x_i, t - \tau_i) = f(s - x_i)}
#' where \eqn{f} is the density function of \eqn{\textrm{N}(0, \Sigma)}.
#' \item For time-dependent spatial excitement (\code{time_independent = FALSE}),
#' \eqn{G_i(s-x_i, t - \tau_i) = f(s - x_i)}
#' where \eqn{f} is the density function of \eqn{\textrm{N}(0, (t-\tau_i)\Sigma)}.
#' }
#' @return A list containing components of the fitted model, see \code{TMB::MakeADFun}.
#' @seealso \code{\link{fit_hawkes}} and \code{\link{fit_lgcp}}
#' @examples \donttest{
#' ## No GMRF
#' if(requireNamespace("INLA")){
#' data(xyt, package = "stelfi")
#' N <- 50
#' locs <- data.frame(x = xyt$x[1:N], y = xyt$y[1:N])
#' times <- xyt$t[1:N]
#' domain <- sf::st_as_sf(xyt$window)
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd, max.edge = 0.75, cutoff = 0.3) 
#' param <- list( mu = 3, alpha = 1, beta = 3, xsigma = 0.2, ysigma = 0.2, rho = 0.8)
#' fit <- fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh, parameters = param) 
#' get_coefs(fit)
#' ## GMRF
#' param <- list( mu = 5, alpha = 1, beta = 3, kappa = 0.9, tau = 1, xsigma = 0.2,
#' ysigma = 0.2, rho = 0.8)
#' fit <- fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh,
#' parameters = param, GMRF = TRUE)
#' get_coefs(fit)
#' }
#' }
#' @export
fit_stelfi <-  function(times, locs, sf, smesh,  parameters, covariates,
                        GMRF = FALSE,
                        time_independent = TRUE,
                        tmb_silent = TRUE,
                        nlminb_silent = TRUE, ...) {
    ## parameters
    coefs <- parameters[["coefs"]]
    if (is.null(coefs)) {
      if (!missing(covariates)) {
        coefs <- numeric(ncol(covariates) + 1)
        coefs[1] <- log(0.5 * length(times) / max(times))
      } else {
        coefs <- c(log(0.5 * length(times) / max(times)))
      }
    }
    logit_abratio <- stats::qlogis(parameters[["alpha"]] / parameters[["beta"]])
    if (is.null(logit_abratio)) {
      logit_abratio <- 0
    }
    log_beta <- log(parameters[["beta"]])
    if (is.null(log_beta)) {
      log_beta <- log(2) - log(coefs[1])
    }
    log_xsigma <-  log(parameters[["xsigma"]])
    if (is.null(log_xsigma)) {
      log_xsigma <- 0
    }
    log_ysigma <-  log(parameters[["ysigma"]])
    if (is.null(log_ysigma)) {
      log_ysigma <- 0
    }
    atanh_rho <- atanh(parameters[["rho"]])
    if (is.null(atanh_rho)) {
      atanh_rho <- 0
    }
    ## error checking
    for (i in 2:length(times)) {
      if ((times[i] - times[i - 1]) < 1.e-10)
        stop("times must be in ascending order with no simultaneous events")
    }
    if (is.na(logit_abratio) || is.null(logit_abratio)) {
      stop("alpha/beta must be between 0 and 1")
    }
    if (is.na(log_beta) || is.null(log_beta)) {
      stop("beta must be positive")
    }
    if(sum(names(locs) %in% c("x","y")) < 2)
      stop("Named variables x and y required in arg locs")
    if (nrow(locs) != length(times)) {
      stop("different number of times and spatial locations")
    }
    ## design matrix and covariates
    if(!missing(covariates)) {
      if(length(coefs) != (ncol(covariates) + 1))
        stop("arg coefs should be length ncol.covariates + 1")
      if(nrow(covariates) != nrow(smesh$loc))
         stop("nrow.covariates should be same as spatial mesh size")
      designmat <- cbind(1, covariates)
    } else {
      if(length(coefs) != 1){
        stop("arg coefs should be length 1 if covariates missing")
      }
      designmat <- matrix(rep(1, smesh$n), ncol = 1)
    }
    ## weights
    w <- get_weights(mesh = smesh, sf = sf, plot = FALSE)$weights
    locs <- as.matrix(locs)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    if(!GMRF) { ## No GMRF
      res <- fit_hspat_tmb(times = times, locs = locs, sf = sf, w  = w,
                            smesh = smesh, coefs = coefs, designmat = designmat,
                            logit_abratio = logit_abratio, log_beta = log_beta,
                            log_xsigma = log_xsigma,
                            log_ysigma = log_ysigma,
                            atanh_rho = atanh_rho,
                            lmat = lmat, simple = as.numeric(time_independent),
                            tmb_silent = tmb_silent,
                            nlminb_silent = nlminb_silent, ...)
    }else{ ## SPDE
        spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
        log_tau <- log(parameters[["tau"]])
        log_kappa <- log(parameters[["kappa"]])
        res <- fit_hspde_tmb(times = times, locs = locs, sf = sf,
                            spde = spde, w  = w,
                            smesh = smesh, coefs = coefs, designmat = designmat,
                            logit_abratio = logit_abratio, log_beta = log_beta,
                            log_kappa = log_kappa,
                            log_tau = log_tau, log_xsigma =  log_xsigma,
                            log_ysigma =  log_ysigma, atanh_rho = atanh_rho,
                            lmat = lmat, simple = as.numeric(time_independent),
                            tmb_silent = tmb_silent,
                            nlminb_silent = nlminb_silent, ...)
    }
    return(res)
}
