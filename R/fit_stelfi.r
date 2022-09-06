#' Function to fit a spde spatial Hawkes process using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
fit_hspde_tmb <- function(times, locs, sp, mesh,
                          coefs, designmat, logit_abratio = 0, log_beta = 0,
                          log_kappa = 0, log_tau = 0, log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0,
                          spde, w, tmax = max(times),
                          reltol = 1e-12, abstol = 1e-12,
                          lmat = lmat, simple = simple,
                          tmb_silent, nlminb_silent, ...) {
    if (!"spde_hawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spde_hawkes")
    }
    sp_sf <- sf::st_as_sf(sp)
    innerloc <- Reduce(rbind, apply(mesh$graph$tv, 1, function(vt) {
        temp <- sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp <- sp::SpatialPolygons(list(temp))
        temp_sf <- sf::st_as_sf(temp)
        if (is.null(sf::st_intersection(temp_sf, sp_sf)))
          rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs,
                xyloc = mesh$loc[, 1:2], reltol = reltol, abstol = abstol,
                spde = spde$param.inla[c("M0", "M1", "M2")], w = w, tmax = tmax,
                designmat = designmat, tv = innerloc, simple = simple, lmat = lmat)
    param <- list(coefs = coefs, logit_abratio = logit_abratio, log_beta = log_beta,
                  log_xsigma =  log_xsigma, log_ysigma =  log_ysigma,
                  atanh_rho = atanh_rho, log_kappa = log_kappa, log_tau = log_tau,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = 1))
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE, random = c("x"),
                           DLL = "spde_hawkes", silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$objective
    return(obj)
}
#' Function to fit a spatial Hawkes process using \code{TMB}
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
#' @inheritParams fit_hspde_tmb
fit_hspat_tmb <- function(times, locs, sp,
                          mesh, coefs, designmat, logit_abratio = 0, log_beta = 0,
                          log_xsigma = 0, log_ysigma = 0, atanh_rho = 0, w,
                          reltol = 1e-12, abstol = 1e-12,
                          tmax = max(times),
                          lmat, simple,
                          tmb_silent, nlminb_silent, ...) {
    if (!"spatial_hawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spatial_hawkes")
    }
    sp_sf <- sf::st_as_sf(sp)
    innerloc <- Reduce(rbind, apply(mesh$graph$tv, 1, function(vt) {
        temp <- sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp <- sp::SpatialPolygons(list(temp))
        temp_sf <- sf::st_as_sf(temp)
        if (is.null(sf::st_intersection(temp_sf, sp_sf)))
          rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs, 
                 xyloc = mesh$loc[,1:2], reltol = reltol, abstol = abstol,
                  w = w, tmax = tmax, designmat = designmat,
                 tv = innerloc, simple = simple, lmat = lmat)
    param <- list(coefs = coefs, logit_abratio = logit_abratio,
                  log_beta = log_beta,
                  log_xsigma =  log_xsigma,
                  log_ysigma =  log_ysigma, atanh_rho = atanh_rho)
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE,
                           DLL = "spatial_hawkes", silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
    obj$objective <- opt$objective
    return(obj)
}
#' Function to fit a spatiotemporal Hawkes model.
#' The self-excitement is Gaussian in space and exponentially decaying in time.
#' @param locs A \code{data.frame} of \code{x} and \code{y} locations, 2xn.
#' @param parameters a list of named parameters
#' Default values used if not provided.
#' "coefs"--logged base rate of the Hawkes process and coefficients of covariates
#' "alpha"--intensity jump after an event occurrence
#' "beta"-- rate of exponential decay of intensity after event occurrence
#' "tau"-- \code{tau} parameter for the GMRF (only if \code{GMRF} == TRUE)
#' "kappa"--\code{kappa} parameter for the GMRF (only if \code{GMRF} == TRUE)
#' @param covariates Optional, a \code{matrix} of covariates at each
#' \code{smesh} node.
#' @param GMRF logical, default FALSE. If TRUE, a spatial Gaussian Markov Random Field is fitted.
#' @param time_independent logical, default TRUE. If FALSE, Gaussian kernels have a
#' covariate matrix that is proportional to time since the event.
#' Warning, this is very memory intensive.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
#' @inheritParams fit_hspat_tmb
#' @export
fit_stelfi <-  function(times, locs, sp, smesh,  parameters, covariates,
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
    w <- get_weights(mesh = smesh, sf = sf::st_as_sf(sp), plot = FALSE)$weights
    locs <- as.matrix(locs)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    if(!GMRF) { ## No GMRF
      res <- fit_hspat_tmb(times = times, locs = locs, sp = sp, w  = w,
                            mesh = smesh, coefs = coefs, designmat = designmat,
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
        res <- fit_hspde_tmb(times = times, locs = locs, sp = sp,
                            spde = spde, w  = w,
                            mesh = smesh, coefs = coefs, designmat = designmat,
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
