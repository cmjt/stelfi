#' Function to fit a spde spatial hawkes process using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
fit_hspde_tmb <- function(times , locs , sp,
                          mesh, lmat,
                          log_mu = 0, logit_abratio = 0, log_beta = 0,
                          log_kappa = 0, log_tau = 0, log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0,
                          spde, w , tmax = max(times),
                          reltol = 1e-12, abstol = 1e-12,
                          ...){
    if (!"spdehawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spdehawkes")
    }
    innerloc <- Reduce(rbind, apply(mesh$graph$tv, 1, function(vt){
        temp = sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp = sp::SpatialPolygons(list(temp))
        if (is.null(rgeos::gIntersection(temp, sp)))
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
                           DLL = "spdehawkes")
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
    return(obj)
}
#' Function to fit a spatial hawkes process using \code{TMB}
#' @inheritParams fit_hawkes
#' @inheritParams fit_lgcp
fit_hspat_tmb <- function(times, locs, sp,
                          mesh,log_mu = 0, logit_abratio = 0, log_beta = 0,
                          log_xsigma = 0,
                          log_ysigma = 0, atanh_rho = 0, w,
                          reltol = 1e-12, abstol = 1e-12,
                          tmax = max(times), ...){
    if (!"spatialhawkes" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("spatialhawkes")
    }
    innerloc = Reduce(rbind, apply(mesh$graph$tv, 1, function(vt){
        temp = sp::Polygons(list(sp::Polygon(mesh$loc[rep(vt, length = length(vt) + 1), 1:2])), '0')
        temp = sp::SpatialPolygons(list(temp))
        if (is.null(rgeos::gIntersection(temp, sp)))
            rep(NULL, 3)
        else
            vt
    }))
    data <- list(times = times, locs = locs, 
                 xyloc = mesh$loc[,1:2], reltol = reltol, abstol = abstol,
                  w = w, tmax = tmax,
                 tv = innerloc)
    param <- list(log_mu = log_mu, logit_abratio = logit_abratio,
                  log_beta = log_beta,
                  log_xsigma =  log_xsigma,
                  log_ysigma =  log_ysigma, atanh_rho = atanh_rho)
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE,
                           DLL = "spatialhawkes")
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
    return(obj)
}
#' Funtion to fit stelfi
fit_stelfi <-  function(times, locs, sp, smesh,  parameters,
                        gaussian = TRUE, ...) {
    ## convert svs
    log_mu <- log(parameters[["mu"]])
    logit_abratio <- stats::qlogis(parameters[["alpha"]] / parameters[["beta"]])
    log_beta <- log(parameters[["beta"]])
    log_xsigma <-  log(parameters[["xsigma"]])
    log_ysigma <-  log(parameters[["ysigma"]])
    atanh_rho <- atanh(parameters[["rho"]])
    ## weights
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)$weights
    if(gaussian == TRUE){
        res <- fit_hspat_tmb(times = times, locs = locs, sp = sp,
                             w  = w,
                             mesh = smesh, log_mu = log_mu,
                             logit_abratio = logit_abratio, log_beta = log_beta,
                             log_xsigma = log_xsigma,
                             log_ysigma = log_ysigma,
                             atanh_rho = atanh_rho, ...)
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
                             log_ysigma =  log_ysigma, atanh_rho = atanh_rho, ...)
    }
    return(res)
}


    
