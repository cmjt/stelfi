#' Simulate a log-Gaussian Cox process (LGCP)
#' 
#' Simulate a realisation of a log-Gaussian Cox process (LGCP) using the
#' \code{TMB} \code{C++} template. If \code{rho} is supplied in \code{parameters}
#' as well as \code{tmesh} then spatiotemporal (AR(1)) data will be simulated.
#'
#' @inheritParams fit_lgcp
#' @param all Logical, if \code{TRUE} then all model components are returned.
#' @return A named list. If \code{all = FALSE} then only the simulated values of
#' the GMRF at each mesh node are returned, \code{x}, alongside the number of
#' events, \code{y}, simulated at each node.
#'
#' @examples \donttest{
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd,
#' max.edge = 0.75, cutoff = 0.3)
#' parameters <- c(beta = 1, log_tau = log(1), log_kappa = log(1))
#' sim <- sim_lgcp(parameters = parameters, sf = domain, smesh = smesh)
#' ## spatiotemporal
#' ndays <- 2
#' w0 <- 2
#' tmesh <- INLA::inla.mesh.1d(seq(0, ndays, by = w0))
#' parameters <- c(beta = 1, log_tau = log(1), log_kappa = log(1), atanh_rho = 0.2)
#' sim <- sim_lgcp(parameters = parameters, sf = domain, smesh = smesh, tmesh = tmesh)
#' }
#' @export
sim_lgcp <- function(parameters, sf,
                     smesh, tmesh,
                     covariates,
                     all = FALSE) {
    ## svs
    beta <- parameters[["beta"]]
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    ## Verify that arguments are correct size and class, and basic processing
    if(!missing(covariates)) {
        if(!"matrix" %in% class(covariates))
            stop("arg covariates must be a matrix")
        if(length(beta) != (ncol(covariates) + 1))
            stop("arg beta should be length ncol.covariates + 1")
    } else {
        if(length(beta) != 1)
            stop("arg beta should be length 1 if covariates missing")
    }
    if (length(log_tau) != 1) stop("log_tau must be a single number")
    if (length(log_kappa) != 1) stop("log_kappa must be a single number")
    ## prepare variables
    if (!missing(tmesh)) {
        atanh_rho <- parameters[["atanh_rho"]]
        if (length(atanh_rho) != 1)
            stop("atanh_rho must be a single number")
        if (!missing(covariates)) if (nrow(covariates) != nrow(smesh$loc) * k)
                                      stop("nrow.covariates should be size of spatial mesh by number of time knots")
        k <- length(tmesh$loc)
        locs <- matrix(0, nrow = k * 10, ncol = 2)
        locs <- data.frame(x = locs[, 1], y = locs[, 2], t = rep(0:(k - 1), each = 10))
        tmp <- prep_data_lgcp(locs = locs, sf = sf, smesh = smesh, tmesh = tmesh)
        
    } else {
        if(!missing(covariates)) {
            if(length(beta) != (ncol(covariates) + 1))
                stop("arg beta should be length ncol.covariates + 1")
        }
        locs <-  matrix(0, nrow = 10, ncol = 2)
        locs <- data.frame(x = locs[, 1], y = locs[, 2])
        k <- 1
        atanh_rho <- NULL
        tmp <- prep_data_lgcp(locs = locs, sf = sf, smesh = smesh)
    }
    ## Designmat
    if(!missing(covariates)) {
        designmat <- cbind(1, covariates)
    } else {
        designmat <- matrix(rep(1, length(tmp$ypp)), ncol = 1)
    }
    ## Model fitting
    res <- fit_lgcp_tmb(y = tmp$ypp, A = tmp$A,
                        designmat = designmat,
                        spde = tmp$spde$param.inla[c("M0", "M1", "M2")],
                        w = tmp$w,
                        idx = tmp$idx, beta = beta,
                        x = matrix(0, nrow = tmp$spde$n.spde, ncol = k),
                        log_tau = log_tau, log_kappa = log_kappa,
                        atanh_rho = atanh_rho, tmb_silent = TRUE,
                        nlminb_silent = TRUE,
                        simulation = TRUE)
    ## Simulation
    res$env$last.par["beta"] <- beta
    res$env$last.par["log_tau"] <- log_tau
    res$env$last.par["log_kappa"] <- log_kappa
    if (!missing(tmesh)) {
        res$env$last.par["atanh_rho"] <- atanh_rho
    }
    simdata <- res$simulate(complete = all)
    return(simdata)
}
