#' Fit a spatial or spatiotemporal log-Gaussian Cox process (LGCP)
#'
#' \code{fit_lgcp_tmb} fits a LGCP using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field. For
#' a simpler wrapper use \code{\link{fit_lgcp}} as this is an internal function.
#' 
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
#' \code{INLA::inla.spde.make.A()}.
#' @param designmat The design matrix for the fixed effects.
#' @param spde The structure of SPDE object as defined in
#' \code{INLA::inla.spde2.matern()}.
#' The minimal required components are \code{M0}, \code{M1}, \code{M2}.
#' @param w A vector of model weights; corresponds to the \code{E} term for
#' poisson models, see \code{INLA::inla.doc("poisson")} for more detail.
#' @param idx A binary vector of the same size as the observation
#' vector \code{\link{y}}. With this vector, the log-likelihood can
#' be computed using a subset of the observations: 1 for contributing
#' to the log-likelihood, and 0 otherwise.
#' @param beta A vector of fixed effects coefficients to be estimated
#' (same length as \code{ncol(designmat)}.
#' @param x The random field/effects. Set this variable random in
#' \code{TMB::MadeADFun}.
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
#' @param atanh_rho optional, \code{arctan(rho)} AR1 parameter.
#' @param simulation \code{logical}, simulate data, default \code{FALSE}.
#' @inheritParams fit_lgcp
#' @noRd
fit_lgcp_tmb <-  function(y, A, designmat, spde, w, idx, beta,
                          log_tau, log_kappa,
                          atanh_rho, x, tmb_silent,
                          nlminb_silent, simulation, ...) {
    data <- list(y = y, A = A, designmat = designmat,
                 spde = spde, w = w,
                 idx = idx,
                 model_type = "lgcp")
    if(is.null(atanh_rho)) {
        param <- list(beta = beta,  log_kappa = log_kappa,
                      log_tau = log_tau,
                      x = x)
    }else{
        param <- list(beta = beta,  log_kappa = log_kappa,
                      log_tau = log_tau,
                      atanh_rho = atanh_rho, x = x)
    }
    obj <- TMB::MakeADFun(data = data, parameters = param, hessian = TRUE,
                          random = c("x"), DLL = "stelfi",
                          silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    if (simulation == FALSE) {
        opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                             control = list(trace = trace), ...)
        obj$objective <- opt$objective
    }
    return(obj)
}
#' Spatial or spatiotemporal log-Gaussian Cox process (LGCP)
#'
#' Fit a log-Gaussian Cox process (LGCP) using Template Model Builder (TMB) and the
#' \code{R_inla} namespace for the SPDE construction of the latent field. 
#'
#' @details A log-Gaussian Cox process (LGCP) where the Gaussian random field, \eqn{Z(\boldsymbol{x})},
#' has zero mean, variance-covariance matrix \eqn{\boldsymbol{Q}^{-1}}, and covariance function
#' \eqn{C_Z}. The random intensity surface is
#' \eqn{\Lambda(\boldsymbol{x}) = \textrm{exp}(\boldsymbol{X}\beta + G(\boldsymbol{x}) + \epsilon)},
#' for design matrix \eqn{\boldsymbol{X}}, coefficients \eqn{\boldsymbol{\beta}}, and random error \eqn{\epsilon}.
#'
#' 
#' Shown in Lindgren et. al., (2011) the stationary solution to the SPDE (stochastic
#' partial differential equation) \eqn{(\kappa^2 - \Delta)^{\frac{\nu + \frac{d}{2}}{2}}G(s) = W(s)} is
#' a random field with a Matérn covariance function,
#' \eqn{C_Z \propto {\kappa || x - y||}^{\nu}K_{\nu}{\kappa || x - y||}}. Here \eqn{\nu} controls
#' the smoothness of the field and \eqn{\kappa} controls the range.
#'
#' 
#' A Markovian random field is obtained when \eqn{\alpha = \nu + \frac{d}{2}} is an integer. Following
#' Lindgren et. al., (2011) we set \eqn{\alpha = 2} in 2D and therefore fix \eqn{\nu = 1}. Under these
#' conditions the solution to the SPDE is a Gaussian Markov Random Field (GMRF). This is the approximation
#' we use.
#'
#' The (approximate) spatial range \eqn{= \frac{\sqrt{8 \nu}}{\kappa} = \frac{\sqrt{8}}{\kappa}} and
#' the standard deviation of the model, \eqn{\sigma = \frac{1}{\sqrt{4 \pi \kappa^2 \tau^2}}}.
#' Under \code{INLA} (Lindgren and Rue, 2015) methodology the practical range is defined as the
#' distance such that the correlation is \eqn{\sim 0.139}.
#'
#'
#' @references Lindgren, F., Rue, H., and Lindström, J. (2011)
#' An explicit link between {G}aussian fields and {G}aussian {M}arkov random fields: the stochastic
#' partial differential equation approach. \emph{Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology)}, \strong{73}: 423--498.
#'
#' @references Lindgren, F. and Rue, H. (2015) Bayesian spatial modelling with R-INLA.
#' \emph{Journal of Statistical Software}, \strong{63}: 1--25.
#' 
#' @param sf An \code{sf} of type \code{POLYGON} specifying the spatial region
#' of the domain.
#' @param locs A \code{data.frame} of \code{x} and \code{y} locations, \eqn{2 \times n}. If a
#' spatiotemporal model is to be fitted then there should be the third column (\code{t}) of the occurance times.
#' @param smesh A Delaunay triangulation of the spatial domain returned by \code{INLA::inla.mesh.2d()}.
#' @param tmesh Optional, a temporal mesh returned by \code{INLA::inla.mesh.1d()}.
#' @param parameters A named list of parameter starting values:
#' \itemize{
#' \item \code{beta}, a vector of fixed effects coefficients to be estimated, \eqn{\beta}
#' (same length as \code{ncol(covariates)} + 1 );
#' \item  \code{log_tau}, the \eqn{\textrm{log}(\tau)} parameter for the GMRF;
#' \item \code{log_kappa}, \eqn{\textrm{log}(\kappa)} parameter for the GMRF;
#' \item \code{atanh_rho}, optional, \eqn{\textrm{arctan}(\rho)} AR1 temporal parameter.
#' }
#' Default values are used if none are provided. NOTE: these may not always be appropriate.
#' @param covariates Optional, a \code{matrix} of covariates at each
#' \code{smesh} and \code{tmesh} node combination.
#' @param tmb_silent Logical, if \code{TRUE} (default) then
#' TMB inner optimisation tracing information will be printed.
#' @param nlminb_silent Logical, if \code{TRUE} (default) then for each iteration
#' \code{nlminb()} output will be printed.
#' @param ... optional extra arguments to pass into \code{stats::nlminb()}.
#' @seealso \code{\link{fit_mlgcp}} and \code{\link{sim_lgcp}}
#' @examples \donttest{
#' ### ********************** ###
#' ## A spatial only LGCP
#' ### ********************** ###
#' data(xyt, package = "stelfi")
#' domain <- sf::st_as_sf(xyt$window)
#' locs <- data.frame(x = xyt$x, y = xyt$y)
#' stelfi_load_inla()
#' bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
#' smesh <- INLA::inla.mesh.2d(boundary = bnd,
#' max.edge = 0.75, cutoff = 0.3)
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
#' parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))
#' ### ********************** ###
#' ## A spatiotemporal LGCP, AR(1)
#' ### ********************** ###
#' ndays <- 2
#' locs <- data.frame(x = xyt$x, y = xyt$y, t = xyt$t)
#' w0 <- 2
#' tmesh <- INLA::inla.mesh.1d(seq(0, ndays, by = w0))
#' fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh, tmesh = tmesh,
#'  parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1), atanh_rho = 0.2))
#' }
#' @export
fit_lgcp <-  function(locs, sf, smesh, tmesh, parameters = list(), covariates,
                      tmb_silent = TRUE,
                      nlminb_silent = TRUE, ...) {
    ## read in parameters
    log_tau <- parameters[["log_tau"]]
    if (is.null(log_tau)) {
      log_tau <- 0
    }
    log_kappa <- parameters[["log_kappa"]]
    if (is.null(log_kappa)) {
      log_kappa <- 0
    }
    beta <- parameters[["beta"]]
    # default value of the intercept is log of the density of points
    if (is.null(beta)) {
      area <- sum(get_weights(smesh, sf)$weights)
      avg_rate <- log(nrow(locs)/area)
      if (!missing(covariates)) {
        beta <- numeric(length= 1 + ncol(covariates))
        beta[1] <- avg_rate
      } else {
        beta <- avg_rate
      }
    }
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
        tmp <- prep_data_lgcp(locs = locs, sf = sf, smesh = smesh, tmesh = tmesh)
        k <- length(tmesh$loc)
        if (!missing(covariates)) if (nrow(covariates) != nrow(smesh$loc)*k)
            stop("nrow.covariates should be size of spatial mesh by number of time knots")
        atanh_rho <- parameters[["atanh_rho"]]
        if (is.null(atanh_rho)) {
          atanh_rho <- 0
        }
        if (length(atanh_rho) != 1)
            stop("atanh_rho must be a single number")
    } else {
        if (!missing(covariates)) if(nrow(covariates) != nrow(smesh$loc))
            stop("nrow.covariates should be same as spatial mesh size")
        tmp <- prep_data_lgcp(locs = locs, sf = sf, smesh = smesh)
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

#' Internal function to prep data as \code{INLA::stack}
#' @inheritParams fit_lgcp
#' @noRd
prep_data_lgcp <- function(locs, sf, smesh, tmesh) {
    ## E
    w <- get_weights(mesh = smesh, sf = sf, plot = FALSE)
    w_areas <- w$weights
    nv <- smesh$n
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## spatial or spatiotemporal
    if (!missing(tmesh)) {
        k <- length(tmesh$loc)
        t.breaks <- sort(c(tmesh$loc[c(1, k)],
                           tmesh$loc[2:k - 1] / 2 + tmesh$loc[2:k] / 2))
        time <- factor(findInterval(locs$t, t.breaks),
                       levels = 1:(length(t.breaks) - 1))
        y.pp <- list()
        for(i in unique(time)) {
            y.pp[[i]] <-  points_in_mesh(locs[time == i, ], w)
        }
        ypp <- as.numeric(unlist(y.pp))
        w.t <- Matrix::diag(INLA::inla.mesh.fem(tmesh)$c0)
        expected <- w_areas[rep(1:nv, k)] * (w.t[rep(1:k, each = nv)])
        A <- INLA::inla.spde.make.A(smesh, smesh$loc[rep(1:nv, k), ],
                              group = rep(1:k, each = nv), mesh.group = tmesh)
        idx <- rep(1, length(ypp))
    }else{
        ypp <- points_in_mesh(locs, w)
        expected <- w_areas
        A <- Matrix::sparseMatrix(i = 1:nv, j = 1:nv, x = 1)
        idx <- expected > 0
    }
    lst <- list(ypp = ypp, A = A, spde = spde, w = expected, idx = idx)
    return(lst)
}
