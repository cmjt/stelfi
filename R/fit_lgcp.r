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
#' @param beta A vector of fixed effects coefficients to be estimated
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
#' @param ... arguments to pass into \code{nlminb()}
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @export
fit_lgcp_tmb <-  function(y, A, designmat, spde, w, idx, beta,
                          log_tau, log_kappa,
                          atanh_rho, x, tmb_silent,
                          nlminb_silent, simulation,...) {
    if (!"lgcp" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("lgcp")
    }
    data <- list(y = y, A = A, designmat = designmat,
                 spde = spde, w = w,
                 idx = idx)
    if(is.null(atanh_rho)){
        param <- list(beta = beta,  log_kappa = log_kappa,
                      log_tau = log_tau,
                      x = x)
    }else{
        param <- list(beta = beta,  log_kappa = log_kappa,
                      log_tau = log_tau,
                      atanh_rho = atanh_rho, x = x)
    }
    obj <- TMB::MakeADFun(data = data, parameters = param,
                          random = c("x"), DLL = "lgcp",
                          silent = tmb_silent)
    obj$hessian <- TRUE
    trace <- if(nlminb_silent) 0 else 1
    if (simulation == FALSE) {
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
        obj$objective <- opt$objective
    }
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
#' @param parameters a list of named parmeters:
#' "beta"--A vector of fixed effects coefficients to be estimated
#' (same length as \code{ncol(covariates)} + 1 )
#'  "log_tau"-- tau parameter for the GMRF
#'  "log_kappa"-- kappa parameter for the GMRF
#'  "rho"-- optional, rho AR1 parameter
#' @param covariates optional, a \code{matrix} of covariates at each
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @export
fit_lgcp <-  function(locs, sp, smesh, tmesh, parameters, covariates,
                      tmb_silent = TRUE,
                      nlminb_silent = TRUE, ...) {
    # Verify that arguments are correct size and class, basic processing
    if(sum(names(locs) %in% c("x","y")) < 2) stop("Named variables x and y required in arg locs")
    if(!missing(covariates)) if(!"matrix" %in% class(covariates)) stop("arg covariates must be a matrix")
    beta <- parameters[["beta"]]
    if(!missing(covariates)) if(length(beta) != (ncol(covariates) + 1))
        stop("arg beta should be length ncol.covariates + 1")
    if (missing(covariates)) if(length(beta) != 1)
        stop("arg beta should be length 1 if covariates missing")
    if (!missing(tmesh)) {
        if(!"t" %in% names(locs)) stop("Need a variable named t in arg locs")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
        k <- length(tmesh$loc)
        if (!missing(covariates)) if (nrow(covariates) != nrow(smesh$loc)*k)
            stop("nrow.covariates should be size of spatial mesh times number of time knots")
        atanh_rho <- parameters[["atanh_rho"]]
    }else{
        if (!missing(covariates)) if(nrow(covariates) != nrow(smesh$loc))
            stop("nrow.covariates should be same as spatial mesh size")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
        k <- 1
        atanh_rho <- NULL
    }
    ## svs
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    ## Designmat
    if(!missing(covariates)) {
        designmat <- cbind(1, covariates)
    } else {
        designmat <- matrix(rep(1, length(tmp$ypp) ), ncol = 1)
    }
    ## Model fitting
    res <- fit_lgcp_tmb(y = tmp$ypp, A = tmp$A,
                        designmat = designmat,
                        spde = tmp$spde$param.inla[c("M0", "M1", "M2")],
                        w = tmp$w,
                        idx = tmp$idx , beta = beta,
                        x = matrix(0, nrow = tmp$spde$n.spde, ncol = k),
                        log_tau = log_tau, log_kappa = log_kappa,
                        atanh_rho = atanh_rho, tmb_silent = tmb_silent,
                        nlminb_silent = nlminb_silent,
                        simulation = FALSE,
                        ...)
    return(res)

}

#' Function to prep data as per INLA stack
#' @inheritParams fit_lgcp
prep_data_lgcp <- function(locs, sp, smesh, tmesh) {
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys
    nv <- smesh$n

    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## spatial or spatiotemporal
    if (!missing(tmesh)) {
        k <- length(tmesh$loc)
        area <- factor(sp::over(sp::SpatialPoints(cbind(xyt$x, xyt$y)), polys),
                       levels = 1:length(polys))

        t.breaks <- sort(c(tmesh$loc[c(1, k)],
                           tmesh$loc[2:k - 1] / 2 + tmesh$loc[2:k] / 2))
        time <- factor(findInterval(xyt$t, t.breaks),
                       levels = 1:(length(t.breaks) - 1))
        agg.dat <- as.data.frame(table(area, time))
        for(j in 1:2) # set time and area as integer
            agg.dat[[j]] <- as.integer(as.character(agg.dat[[j]]))
        ypp <- agg.dat$Freq
        w.t <- Matrix::diag(INLA::inla.mesh.fem(tmesh)$c0)
        expected <- w_areas[agg.dat$area] * (w.t[agg.dat$time])
        A <- INLA::inla.spde.make.A(smesh, smesh$loc[agg.dat$area, ],
                              group = agg.dat$time, mesh.group = tmesh)
        idx <- rep(1, length(ypp))
    }else{
        ypp <- stelfi::points.in.mesh(locs, polys)
        expected <- w_areas
        A <- Matrix::sparseMatrix(i = 1:nv, j = 1:nv, x = 1)
        idx <- expected > 0
    }
    lst <- list(ypp = ypp, A = A, spde = spde, w = expected, idx = idx)
    return(lst)
}
# FUNCTION HEADER TO BE WRITTEN
#' Function to fit a spatial or spatiotemporal log-Gaussian Cox process using \code{TMB}
#'
#' A simple to use wrapper for \code{fit_lgcp_tmb}
#' @param sp \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} of the domain
#' @param locs 2xn \code{data.frame} of locations x, y. If locations have
#' time stapms then this is the third column of the 3xn matrix
#' @param smesh spatial mesh of class \code{"inla.mesh"}
#' @param tmesh optional, temporal mesh of class \code{"inla.mesh.1d"}
#' @param parameters a list of named parmeters:
#' "beta"--A vector of fixed effects coefficients to be estimated
#' (same length as \code{ncol(covariates)} + 1 )
#'  "log_tau"-- tau parameter for the GMRF
#'  "log_kappa"-- kappa parameter for the GMRF
#'  "rho"-- optional, rho AR1 parameter
#' @param covariates optional, a \code{matrix} of covariates at each
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @export
simulate_lgcp <- function(parameters, sp, smesh, covariates) {
    # Currently written to simulate spatial models only
    if(!missing(covariates)) if(!"matrix" %in% class(covariates)) stop("arg covariates must be a matrix")
    if(!missing(covariates)) if(!"matrix" %in% class(covariates)) stop("arg covariates must be a matrix")
    beta <- parameters[["beta"]]
    if(!missing(covariates)) if(length(beta) != (ncol(covariates) + 1))
        stop("arg beta should be length ncol.covariates + 1")
    
    # prepare variables
    locs = matrix(0, nrow=10, ncol=2)
    tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
    k <- 1
    atanh_rho <- NULL
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    ## Designmat
    if(!missing(covariates)) {
        designmat <- cbind(1, covariates)
    } else {
        designmat <- matrix(rep(1, length(tmp$ypp) ), ncol = 1)
    }
    ## Model fitting
    res <- fit_lgcp_tmb(y = tmp$ypp, A = tmp$A,
                        designmat = designmat,
                        spde = tmp$spde$param.inla[c("M0", "M1", "M2")],
                        w = tmp$w,
                        idx = tmp$idx , beta = beta,
                        x = matrix(0, nrow = tmp$spde$n.spde, ncol = k),
                        log_tau = log_tau, log_kappa = log_kappa,
                        atanh_rho = atanh_rho, tmb_silent = TRUE,
                        nlminb_silent = TRUE,
                        simulation = TRUE)
                        #...)
    # Perform simulation
    length_beta <- length(beta)
    length_x <- tmp$spde$n.spde
    res$env$last.par[1:length_beta] <- beta
    res$env$last.par[(length_beta+1):(length_beta+length_x)] <- x
    res$env$last.par[(length_beta+length_x+1)] <- log_tau
    res$env$last.par[(length_beta+length_x+2)] <- log_kappa
    simdata <- res$simulate(complete=TRUE)
    return(simdata)
}
    
    
