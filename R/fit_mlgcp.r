#' Function to fit a spatial or spatiotemporal marked log-Gaussian Cox
#' process using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' For a simpler wrapper use \code{fit_mlgcp()}
#' @param ypp  A vector for point process response (poisson) (see \code{fit_lgcp()})
#' @param marks  A matrix of marks for each observation of the point pattern
#' @param lmat A sparse matrix mapping mesh points to the observations
#' @param spde The structure of spde object as defined in \code{r-inla}.
#' The minimal required components are M0, M1, M2.
#' @param w a vector of model weights; corresponds to the E term for
#' poisson models in \code{INLA},
#' see \code{INLA::inla.doc("poisson")} for more detail.
#' @param idx  a matrix of size (no of resp) * (no of resp + 1).
#' The first column is for sharing random field with point process.
#' For the first column, the values are indicated as >0 if sharing random field with point process.
#' For the i-th column (idx.col(i + 1)), the values are also indicated as >0 if sharing random
#' field defined in i-th random field. If idx(i, i + 1) = 0, then the entire column is ignored.
#' @param strfixed a matrix of fixed structural parameters. For some distributions,
#' non-NaN values will overwrites the corresponding strparam (consider it fixed),
#' and NaN values indicate that it should be estimated (there is no NA in C++).
#' Non-NaN values can be set partially. The entire column shares one estimate if
#' it needs estimating (and hence strparam is a vector of length no. of resp).
#' Details:  Normal distribution, this is the log of element-wise standard deviation.
#' If non-NaN, it is taken as a fixed value  of the log of standard deviation.
#' Poisson distribution, this is the effort. The default input passed from R should be 1s.
#' Binomial distribution, this is the number of trials.
#' The default input from R should be 1s, i.e., the bernoulli distribution.
#' Gamma distribution, this is the log of scale. If non-NaN, this is taken as
#' a fixed value of the log of the scale.
#' @param methods integer:
#' 0 - normal distribution, strparam/strfixed as log_sigma.
#' 1 - poisson distribution, strparam is not referenced, strfixed as effort.
#' 2 - binomial distribution, strparam is not referenced, strfixed as number of trials.
#' 3 - gamma distribution, the implementation in TMB is shape-scale. strparam/strfixed as log_scale;
#' @param betaresp starting value of the intercept term of each response
#' @param betapp starting value of the  intercept of point process
#' @param beta A matrix of the same size as idx.
#' This matrix contains the estimates for the sharing effects.
#' Note that for the right no. of resp-by-no. of resp matrix, the diagonal terms WILL always be 0.
#' This might be resulting from the features of TMB. Consider a manual fixed from R side if needed.
#' All the values in this matrix is constrained to be a normal variable.
#' @param log_kappa  log of kappas for the random field. The lengths of these vectors are no. of resp + 1.
#' The first element is for the random field of the point process.
#' @param log_tau  log of taus for the random field. The lengths of these vectors are no. of resp + 1.
#' The first element is for the random field of the point process.
#' @param strparam see \code{strfixed} and \code{methods}
#' @param silent logical, by default FALSE. If TRUE model fitting progress
#' not printed to console.
#' @param ... arguments to pass into \code{nlminb()}
#' @export
fit_mlgcp_tmb <- function(ypp, marks, lmat, spde, w, idx, strfixed, methods,
                          betaresp, betapp, beta, log_kappa, log_tau,
                          strparam, ...){
    if (!"prefsampling" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("prefsampling")
    }
    data <- list(yresp = marks, ypp = ypp, lmat = lmat,
                 spde = spde$param.inla[c("M0", "M1", "M2")], w = w,
                 idx = idx, methods = methods,
                 strfixed = strfixed)
    param <- list(betaresp = betaresp, betapp = betapp, beta = beta,
                  log_kappa = log_kappa, log_tau = log_tau, strparam = strparam,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = sum(diag(idx[, -1]) > 0) + 1))
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE,
                           random = c("x"), DLL = "prefsampling")
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
    return(obj)
}
#' Function to fit a spatial  marked
#' log-Gaussian Cox process using \code{TMB}
#'
#' A simple to use wrapper for \code{fit_mlgcp_tmb}
#' @param sp \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} of the domain
#' @param locs 2xn \code{data.frame} of locations x, y. If locations have
#' time stapms then this is the third column of the 3xn matrix
#' @param smesh spatial mesh of class \code{"inla.mesh"}
#' @param parameters a list of named parmeters:
#' "beta"--A vector of  coefficients to be estimated
#' of length \code{ncol(marks)} + 1.
#'  "tau"-- tau parameter for the GMRF.
#'  "kappa"-- kappa parameter for the GMRF.
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @export
fit_mlgcp <-  function(locs, sp, marks, smesh, parameters, methods,
                       strfixed, strparam, idx, ...) {
    ## convert svs
    beta <- parameters[["beta"]]
    log_tau <- log(parameters[["tau"]])
    log_kappa <- log(parameters[["kappa"]])
    betaresp <- parameters[["betaresp"]]
    betapp <- parameters[["betapp"]]
    ## data
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys
    points.in.mesh = function(xy, dmesh){
        sapply(1:length(dmesh), function(i){
            coord = raster::geom(dmesh[i, ])[,c("x", "y")]
            sum(sp::point.in.polygon(xy[, 1], xy[, 2], coord[, 1], coord[, 2]) > 0)
        })
    }
    ypp <- points.in.mesh(locs, polys)
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    ## Model fitting
    res <- fit_mlgcp_tmb(ypp = ypp, marks = marks, lmat = lmat,
                         spde = spde,
                         w = w_areas, idx = idx, strfixed = strfixed,
                         methods = methods,
                         betaresp = betaresp, betapp = betapp,
                         beta = beta, log_kappa = log_kappa, log_tau = log_tau,
                         strparam =  strparam,  ...)
    return(res)
}

