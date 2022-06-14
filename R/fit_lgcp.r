#' Fit a spatial or spatiotemporal log-Gaussian Cox process (LGCP)
#'
#' \code{fit_lgcp_tmb} fits a LGCP using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field. For
#' a simpler wrapper use \code{\link{fit_lgcp}} as this is an internal function.
#' @param y The vector of observations.
#' For spatial only, this is the vector of counts.
#' For spatial + AR1 temporal, this vector needs to be
#' arranged in the following way:
#'\tabular{rr}{ y \tab  t \cr
#' y_1  \tab             1 \cr
#' y_2  \tab             1 \cr
#' ... \tab ... \cr
#' y_n  \tab             1 \cr
#' y_{n+1}  \tab         2 \cr
#' y_{n+2}  \tab         2 \cr
#' ... \tab ... \cr
#' y_{n * t_n - 1} \tab  t_n \cr
#' y_{n * t_n}  \tab     t_n \cr
#' }
#' @param A The predictor matrix A, obtained from
#' \code{\link[INLA]{inla.spde.make.A}}.
#' @param designmat The design matrix for the fixed effects.
#' @param spde The structure of SPDE object as defined in
#' \code{\link[INLA]{inla.spde2.matern}}.
#' The minimal required components are \code{M0}, \code{M1}, \code{M2}.
#' @param w A vector of model weights; corresponds to the \code{E} term for
#' poisson models, see \code{\link[INLA]{inla.doc("poisson")}} for more detail.
#' @param idx A binary vector of the same size as the observation
#' vector \code{\link{y}}. With this vector, the log-likelihood can
#' be computed using a subset of the observations: 1 for contributing
#' to the log-likelihood, and 0 otherwise.
#' @param beta A vector of fixed effects coefficients to be estimated
#' (same length as \code{ncol(\link{designmat})}.
#' @param x The random field/effects. Set this variable random in
#' \code{\link[TMB]{MadeADFun}}.
#' The array is of size \code{n} of random effect for each time knot
#' (x.rows()) by number of temporal
#' knots (x.cols()). Therefore, the total number of random effects is
#' x.rows() * x.cols().
#' Each column represents a time knot:
#' \tabular{rrrrrr}{
#'  t     \tab              1    \tab       2   \tab        3   \tab        ...    \tab       t_n \cr
#' random effects \tab x.col(0)  \tab   x.col(1) \tab   x.col(2)  \tab     ...   \tab      x.col(t_n - 1)
#' }
#' @param log_tau \code{log(tau)} parameter for the GMRF.
#' @param log_kappa \code{log(kappa)} parameter for the GMRF.
#' @param atanh_rho Optional, \code{arctan(rho)} AR1 parameter.
#' @param ... Optional extra arguments to pass into \code{\link[stats]{nlminb}}.
#' @param tmb_silent Logical, default \code{TRUE}:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent Logical, default \code{TRUE}:
#' print function and parameters every iteration.
#' @seealso \code{\link{fit_lgcp}}
fit_lgcp_tmb <-  function(y, A, designmat, spde, w, idx, beta,
                          log_tau, log_kappa,
                          atanh_rho, x, tmb_silent,
                          nlminb_silent, simulation, ...) {
    if (!"lgcp" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("lgcp")
    }
    data <- list(y = y, A = A, designmat = designmat,
                 spde = spde, w = w,
                 idx = idx)
    if(is.null(atanh_rho)) {
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
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                             control = list(trace = trace), ...)
        obj$objective <- opt$objective
    }
    return(obj)
}
#' Fit a spatial or spatiotemporal log-Gaussian Cox process (LGCP)
#' 
#' \code{fit_lgcp} fits a LGCP using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field. Ths is
#' the user friendly wrapper for the internal function \code{\link{fit_lgcp_tmb}}. 
#' @seealso \code{\link{fit_lgcp_tmb}}.
#' @param sp A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' of the domain.
#' @param locs A \code{data.frame} of \code{x} and \code{y} locations, 2xn. If locations have
#' time stamps then this should be the third column of the supplied 3xn matrix.
#' @param smesh A spatial mesh of class \code{\link[INLA]{inla.mesh.2d}}.
#' @param tmesh Optional, a temporal mesh of class \code{\link[INLA]{inla.mesh.1d}}.
#' @param parameters A named list of parameters:
#' \code{beta}, a vector of fixed effects coefficients to be estimated
#' (same length as \code{ncol(covariates)} + 1 );
#' \code{log_tau}, the \code{log(tau)} parameter for the GMRF;
#' \code{log_kappa}, \code{log(kappa)} parameter for the GMRF;
#' \code{atanh_rho}, optional, \code{arctan(rho)} AR1 temporal parameter.
#' @param covariates Optional, a \code{matrix} of covariates at each
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @return A fitted \code{\link[TMB]{MakeADFun}} object.
#' @examples \dontrun{
#' data(xyt, package = "stelfi")
#' domain <- as(xyt$window, "SpatialPolygons")
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' smesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain),
#' max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sp = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' }
#' \dontrun{
#' ndays <- 2
#' locs <- data.frame(x = xyt$x, y = xyt$y, t = xyt$t)
#' w0 <- 2
#' tmesh <- INLA::inla.mesh.1d(seq(0, ndays, by = w0))
#' fit <- fit_lgcp(locs = locs, sp = domain, smesh = smesh, tmesh = tmesh,
#'  parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1), atanh_rho = 0.2))
#' 
#' }
#' @export
fit_lgcp <-  function(locs, sp, smesh, tmesh, parameters, covariates,
                      tmb_silent = TRUE,
                      nlminb_silent = TRUE, ...) {
    ## svs
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    beta <- parameters[["beta"]]
    # Verify that arguments are correct size and class, and basic processing
    if(sum(names(locs) %in% c("x","y")) < 2)
        stop("Named variables x and y required in arg locs")
    if(!missing(covariates)) {
        if(!"matrix" %in% class(covariates))
            stop("arg covariates must be a matrix")
        if(length(beta) != (ncol(covariates) + 1))
            stop("arg beta should be length ncol.covariates + 1")
    } else {
        if(length(beta) != 1)
            stop("arg beta should be length 1 if covariates missing")
    }
    if (length(log_tau) != 1)
        stop("log_tau must be a single number")
    if (length(log_kappa) != 1)
        stop("log_kappa must be a single number")
    if (!missing(tmesh)) {
        if(!"t" %in% names(locs))
            stop("Need a variable named t in arg locs")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
        k <- length(tmesh$loc)
        if (!missing(covariates)) if (nrow(covariates) != nrow(smesh$loc)*k)
            stop("nrow.covariates should be size of spatial mesh by number of time knots")
        atanh_rho <- parameters[["atanh_rho"]]
        if (length(atanh_rho) != 1)
            stop("atanh_rho must be a single number")
    } else {
        if (!missing(covariates)) if(nrow(covariates) != nrow(smesh$loc))
            stop("nrow.covariates should be same as spatial mesh size")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
        k <- 1
        atanh_rho <- NULL
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
                        atanh_rho = atanh_rho, tmb_silent = tmb_silent,
                        nlminb_silent = nlminb_silent,
                        simulation = FALSE,
                        ...)
    return(res)

}

#' Internal function to prep data as per \code{link{INLA}} stack
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
        ypp <- points.in.mesh(locs, polys)
        expected <- w_areas
        A <- Matrix::sparseMatrix(i = 1:nv, j = 1:nv, x = 1)
        idx <- expected > 0
    }
    lst <- list(ypp = ypp, A = A, spde = spde, w = expected, idx = idx)
    return(lst)
}
#' Simulate a log-Gaussian Cox process (LGCP)
#' 
#' \code{simulate_lgcp} simulates a LGCP using the \code{TMB} \code{C++}
#' template. If \code{rho} is supplied in \code{parameters}
#' as well as \code{tmesh} then times knots will also be returned.
#' @inheritParams fit_lgcp
#' @param all Logical, if \code{TRUE} then all model components are returned.
#' @return A named list. If \code{all = FALSE} then only the simulated values of 
#' the GMRF at each mesh node are returned, \code{x}, alongside the number of 
#' events, \code{y}, simulated at each node.
#'
#' @examples \dontrun{
#' data(xyt, package = "stelfi")
#' domain <- as(xyt$window, "SpatialPolygons")
#' smesh <- INLA::inla.mesh.2d(boundary = INLA::inla.sp2segment(domain), 
#' max.edge = 0.75, cutoff = 0.3)
#' parameters <- c(beta = 1, log_tau = log(1), log_kappa = log(1))
#' simulate_lgcp(parameters = parameters, sp = domain, smesh = smesh)
#' }
#' \dontrun{
#' ndays <- 2
#' w0 <- 2
#' tmesh <- INLA::inla.mesh.1d(seq(0, ndays, by = w0))
#' parameters <- c(beta = 1, log_tau = log(1), log_kappa = log(1), atanh_rho = 0.2)
#' simulate_lgcp(parameters = parameters, sp = domain, smesh = smesh, tmesh = tmesh)
#' }
#' @export
simulate_lgcp <- function(parameters, sp, smesh, tmesh, covariates,
                          all = FALSE) {
    ## svs
    beta <- parameters[["beta"]]
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    ## Verify that arguments are correct size and class, and basic processing
    if(sum(names(locs) %in% c("x", "y")) < 2)
        stop("Named variables x and y required in arg locs")
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
        locs <- matrix(0, nrow = 10, ncol = 3)
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
        atanh_rho <- parameters[["atanh_rho"]]
        k <- length(tmesh$loc)
        if (!missing(covariates)) if (nrow(covariates) != nrow(smesh$loc) * k)
                                      stop("nrow.covariates should be size of spatial mesh by number of time knots")
        if (length(atanh_rho) != 1)
            stop("atanh_rho must be a single number")
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
    } else {
        locs <-  matrix(0, nrow = 10, ncol = 2)
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
        k <- 1
        if(!missing(covariates)) {
            if(length(beta) != (ncol(covariates) + 1))
                stop("arg beta should be length ncol.covariates + 1")
        }
        k <- 1
        atanh_rho <- NULL
        tmp <- prep_data_lgcp(locs = locs, sp = sp, smesh = smesh)
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
