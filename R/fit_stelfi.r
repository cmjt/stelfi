#' Function to fit a spde spatial hawkes process using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
fit_hspde_tmb <- function(times, locs, sp,
                          mesh, lmat,
                          log_mu = 0, logit_abratio = 0, log_beta = 0,
                          log_kappa = 0, log_tau = 0, log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0,
                          spde, w , tmax = max(times),
                          reltol = 1e-12, abstol = 1e-12,
                          tmb_silent,
                          nlminb_silent, ...){
    if (!"spde_hawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spde_hawkes")
    }
    sp_sf = sf::st_as_sf(sp)
    innerloc <- Reduce(rbind, apply(mesh$graph$tv, 1, function(vt){
        temp = sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp = sp::SpatialPolygons(list(temp))
        temp_sf = sf::st_as_sf(temp)
        if (is.null(sf::st_intersection(temp_sf, sp_sf)))
          rep(NULL, 3)
        else
            vt
    }))
    lmat <- INLA::inla.spde.make.A(mesh, locs)
    data <- list(times = times, locs = locs, 
                xyloc = mesh$loc[,1:2], reltol = reltol, abstol = abstol,
                spde = spde$param.inla[c("M0", "M1", "M2")], w = w, tmax = tmax,
                tv = innerloc, lmat = lmat)
    param <- list(log_mu = log_mu, logit_abratio = logit_abratio, log_beta = log_beta,
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
#' Function to fit a spatial hawkes process using \code{TMB}
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
fit_hspat_tmb <- function(times, locs, sp,
                          mesh,log_mu = 0, logit_abratio = 0, log_beta = 0,
                          log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0, w,
                          reltol = 1e-12, abstol = 1e-12,
                          tmax = max(times),
                          simple,
                          tmb_silent, nlminb_silent, ...){
    if (!"spatial_hawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spatial_hawkes")
    }
    sp_sf = sf::st_as_sf(sp)
    innerloc = Reduce(rbind, apply(mesh$graph$tv, 1, function(vt){
        temp = sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp = sp::SpatialPolygons(list(temp))
        temp_sf = sf::st_as_sf(temp)
        if (is.null(sf::st_intersection(temp_sf, sp_sf)))
          rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs, 
                 xyloc = mesh$loc[,1:2], reltol = reltol, abstol = abstol,
                  w = w, tmax = tmax,
                 tv = innerloc, simple = simple)
    param <- list(log_mu = log_mu, logit_abratio = logit_abratio,
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
#' Funtion to fit stelfi
#' @inheritParams fit_hspat_tmb
#' @inheritParams fit_hspde_tmb
fit_stelfi <-  function(times, locs, sp, smesh,  parameters,
                        gaussian = TRUE,
                        simple = 0, 
                        tmb_silent = TRUE,
                        nlminb_silent = TRUE, ...) {
    ## parameters
    log_mu <- log(parameters[["mu"]])
    if (is.null(log_mu)) {
      log_mu <- log(0.5 * length(times)/max(times))
    }
    logit_abratio <- stats::qlogis(parameters[["alpha"]] / parameters[["beta"]])
    if (is.null(logit_abratio)) {
      logit_abratio <- 0
    }
    log_beta <- log(parameters[["beta"]])
    if (is.null(log_beta)) {
      log_beta <- log(2) - log_mu
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
    
    if(!"matrix" %in% class(locs)) {
      stop("arg locs must be a matrix")
    }
    
    if (nrow(locs) != length(times)) {
      stop("different number of times and spatial locations")
    }
    
    
    ## weights
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)$weights
    if(gaussian == TRUE){
        res <- fit_hspat_tmb(times = times, locs = locs, sp = sp,
                             w  = w,
                             mesh = smesh, log_mu = log_mu,
                             logit_abratio = logit_abratio, log_beta = log_beta,
                             log_xsigma = log_xsigma,
                             log_ysigma = log_ysigma,
                             atanh_rho = atanh_rho,
                             simple = simple,
                             tmb_silent = tmb_silent,
                             nlminb_silent = nlminb_silent, ...)
    }else{
        ## SPDE
        spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
        log_tau <- log(parameters[["tau"]])
        log_kappa <- log(parameters[["kappa"]])
          res <- fit_hspde_tmb(times = times, locs = locs, sp = sp,
                             spde = spde, w  = w,
                             mesh = smesh, log_mu = log_mu,
                             logit_abratio = logit_abratio, log_beta = log_beta,
                             log_kappa = log_kappa,
                             log_tau = log_tau, log_xsigma =  log_xsigma,
                             log_ysigma =  log_ysigma, atanh_rho = atanh_rho,
                             tmb_silent = tmb_silent,
                             nlminb_silent = nlminb_silent, ...)
    }
    return(res)
}
