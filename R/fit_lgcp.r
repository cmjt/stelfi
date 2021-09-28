#' Function to fit a spatial or spatiotemporal log-Gaussian Cox process using
#' \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field. For
#' a simpler wrapper use \code{fit_lgcp()}
#' @param y the vector of observations.
#' For spatial only, this is the vector of counts.
#' For spatial + AR1 temporal, this vector needs to be
#' arranged in the following way:
#' y                 t
#' ---------------------
#' y_1               1
#' y_2               1
#' ...
#' y_n               1
#' y_{n+1}           2
#' y_{n+2}           2
#' ...
#' y_{n * t_n - 1}   t_n
#' y_{n * t_n}       t_n
#' @param A the predicator matrix A, obtained from \code{r-inla}
#' @param designmat the design matrix for the fixed effects
#' @param spde the structure of spde object as defined in \code{r-inla}.
#' The minimal required components are M0, M1, M2
#' @param w a vector of model weights; corresponds to the E term for
#' poisson models in \code{INLA},
#' see \code{INLA::inla.doc("poisson")} for more detail.
#' @param idx a binary vector of the same size as the observation
#' vector \code{y}.
#' With this vector, the log-likelihood can be computed using a subset
#' of the observations. 1 for contributing to the log-likelihood,
#' and 0 otherwise.
#' @param beta A vector of fixed effects needs to be estimated
#' (same length as \code{ncol(designmat)}
#' @param x the random field/effects
#' This is the random field/effects. Set this variable random in
#' \code{MadeADFun()}.
#' This array is of size no. of random effect for each time knots
#' (x.rows()) by no. of temporal
#' knots (x.cols()), and hence the total number of random effects is
#' x.rows() * x.cols().
#' Each column represents a time knot:
#' t
#'                     1           2           3           ...           t_n
#' ----------------------------------------------------------------------
#' random
#' effects x.col(0)     x.col(1)    x.col(2)       ...         x.col(t_n - 1)
#' ----------------------------------------------------------------------
#' @param log_tau \code{log(tau)} parameter for the GMRF
#' @param log_kappa \code{log(kappa)} parameter for the GMRF
#' @param atanh_rho optional, \code{arctan(rho)} AR1 parameter
#' @param silent logical, by default FALSE. If TRUE model fitting progress
#' not printed to console.
#' @param ... arguments to pass into \code{nlminb()}
#' @export
fit_lgcp_tmb <-  function(y, A, designmat, spde, w, idx, beta, x, log_tau, log_kappa,
                          atanh_rho, silent = FALSE,...) {
    if (!"lgcp" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("lgcp")
    }
    data <- list(y = y, A = A, designmat = designmat,
                 spde = spde, w = w,
                 idx = idx)
    param <- list(beta = beta, x = x, log_tau = log_tau,
                  log_kappa = log_kappa, atanh_rho = atanh_rho)
    obj <- TMB::MakeADFun(data = data, parameters = param,
                          random = c("x"), DLL = "lgcp",
                          silent = silent)
    obj$hessian <- TRUE
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
    return(obj)
}
#' Function to fit a spatial or spatiotemporal log-Gaussian Cox process using \code{TMB}
#'
#' A simple to use wrapper for \code{fit_lgcp_tmb}
#' @param sp \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} of the domain
#' @param locs 2xn \code{data.frame} of locations x, y. If locations have
#' time stapms then this is the third column of the 3xn matrix
#' @param smesh spatial mesh of class \code{"inla.mesh"}
#' @param tmesh optional, temporal mesh of class \code{"inla.mesh.1d"}
#' @param tau tau parameter for the GMRF
#' @param kappa kappa parameter for the GMRF
#' @param rho optional, rho AR1 parameter
#' @param covariates optional, a \code{matrix} of covariates at each
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @export
fit_lgcp <-  function(locs, sp, smesh, tmesh, beta, tau, kappa, rho, covariates, silent, ...) {
    if(missing(silent)) silent <- FALSE
    if(sum(names(locs) %in% c("x","y")) < 2) stop("Named variables x and y required in arg locs")
    if(!missing(covariates)) if(!"matrix" %in% class(covariates)) stop("arg covariates must be a matrix")
    if(!missing(covariates)) if(length(beta) != (ncol(covariates) + 1))
        stop("arg beta should be length ncol.covariates + 1")
    if (!missing(tmesh)) {
        if(!"t" %in% names(locs)) stop("Need a variable named t in arg locs")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
        k <- length(tmesh$loc)
        atanh_rho <- atanh(rho)
    }else{
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
        k <- 1
        atanh_rho <- 0
    }
    ## convert svs
    log_tau <- log(tau)
    log_kappa <- log(kappa)
    ## SPDE
    stk <- tmp[[1]]
    a_st <- tmp[[2]]
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## Designmat
    designmat <- matrix(1, nrow = length(stk$data$data$y), ncol = 1)
    if(!missing(covariates)) designmat <- cbind(1, covariates)
    ## Model fitting
    res <- fit_lgcp_tmb(y = stk$data$data$y, A = a_st,
                        designmat = designmat,
                        spde = spde$param.inla[c("M0", "M1", "M2")],
                        w = stk$data$data$exposure,
                        idx = rep(1, length(stk$data$data$y)), beta = beta,
                        x = matrix(0, nrow = spde$n.spde, ncol = k),
                        log_tau = log_tau, log_kappa = log_kappa,
                        atanh_rho = atanh_rho, silent = silent, ...)
    return(res)

}

#' Function to prep data as per INLA stack
#' @inheritParams fit_lgcp
prep_data_lgcp <- function(locs, sp, smesh, tmesh) {
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys
    area <- factor(sp::over(sp::SpatialPoints(cbind(locs$x, locs$y)), polys),
                   levels = seq(1, length(polys)))
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## spatial or spatiotemporal
    if (!missing(tmesh)) {
        k <- length(tmesh$loc)
        idx <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = k)
        w_t <- INLA::inla.mesh.fem(tmesh)$c0
        t_breaks <- sort(c(tmesh$loc[c(1, k)],
                           tmesh$loc[2:k - 1] / 2 + tmesh$loc[2:k] / 2))
        time <- factor(findInterval(locs$t, t_breaks),
                       levels = 1:(length(t_breaks) - 1))
        agg_dat <- as.data.frame(table(area, time))
        for (j in 1:2) # set time and area as integer
            agg_dat[[j]] <- as.integer(as.character(agg_dat[[j]]))
        e0 <- w_areas[agg_dat$area] * (w_t[agg_dat$time])
        a_st <- INLA::inla.spde.make.A(smesh, smesh$loc[agg_dat$area, ],
                         group = agg_dat$time, mesh.group = tmesh)
    }else{
        idx <- INLA::inla.spde.make.index("s", spde$n.spde)
        agg_dat <- as.data.frame(table(area))
        agg_dat[[1]] <- as.integer(as.character(agg_dat[[1]]))
        e0 <- w_areas[agg_dat$area]
        a_st <- INLA::inla.spde.make.A(smesh, smesh$loc[agg_dat$area, ])
    }
    stk <- INLA::inla.stack(
                     data = list(y = agg_dat$Freq, exposure = e0),
                     A = list(1, a_st),
                     effects = list(idx,b0 = rep(1, nrow(agg_dat))))
    
    
    return(list(stk,a_st))
}
