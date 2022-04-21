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
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{nlminb()}
#' @export
fit_mlgcp_tmb <- function(ypp, marks, lmat, spde, w, idx, strfixed, methods,
                          betaresp, betapp, beta_coefs_pp, beta, log_kappa, log_tau,
                          strparam, fields, tmb_silent,
                          nlminb_silent, ...){
    if (!"prefsampling" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("prefsampling")
    }
    data <- list(yresp = marks, ypp = ypp, lmat = lmat,
                 spde = spde$param.inla[c("M0", "M1", "M2")], w = w,
                 idx = idx, methods = methods,
                 strfixed = strfixed, mark_field = fields)
    param <- list(betaresp = betaresp, betapp = betapp, beta = beta, log_kappa = log_kappa, 
                  log_tau = log_tau, strparam = strparam, beta_coefs_pp = beta_coefs_pp,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = sum(diag(idx[, -1]) > 0) + 1))
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE,
                           random = c("x"), DLL = "prefsampling",
                           silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(trace = trace), ...)
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
#' @param fields whether each mark has it's own GMRF
#' @inheritParams fit_lgcp_tmb
#' @export
fit_mlgcp <-  function(locs, sp, marks, smesh, parameters, methods,
                       strfixed, strparam, idx, fields,
                       tmb_silent = TRUE,
                       nlminb_silent = TRUE, ...) {
    
    ## convert svs
    beta <- parameters[["beta"]]
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    betaresp <- parameters[["betaresp"]]
    betapp <- parameters[["betapp"]]
    beta_coefs_pp <- parameters[["beta_coefs_pp"]]
    
    ## Verify args are correct size and class
    if(!"matrix" %in% class(idx)) stop("arg idx must be a matrix")
    n_marks <- ncol(marks)
    if (nrow(idx) != n_marks) stop("nrow.idx must equal ncol.marks")
    if (ncol(idx) != (n_marks+1)) stop("ncol.idx must equal ncol.marks+1")
    if (nrow(beta) != n_marks) stop("nrow.beta must equal ncol.marks")
    if (ncol(beta) != (n_marks+1)) stop("ncol.beta must equal ncol.marks+1")
    if (length(log_tau) != (n_marks+1)) stop("log_tau must have length ncol.marks+1")
    if (length(log_kappa) != (n_marks+1)) stop("log_kappa must have length ncol.marks+1")
    if (length(betaresp) != (n_marks)) stop("betaresp must have length ncol.marks")
    if (length(beta_coefs_pp) != (n_marks)) stop("beta_coefs_pp must have length ncol.marks")
    if (length(betapp) != 1) stop("betapp must have length 1")
    if (length(methods) != n_marks) stop("arg methods must have length ncol.marks")
    if (length(fields) != n_marks) stop("arg fields must have length ncol.marks")
    if (sum(fields) != sum(diag(idx[, -1]) > 0)) stop("idx and/or fields incorrectly set")
    ## data
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys

    ypp <- stelfi::points.in.mesh(locs, polys)
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    ## Model fitting
    res <- fit_mlgcp_tmb(ypp = ypp, marks = marks, lmat = lmat,
                         spde = spde,
                         w = w_areas, idx = idx, strfixed = strfixed,
                         methods = methods,
                         betaresp = betaresp, betapp = betapp, beta_coefs_pp = beta_coefs_pp,
                         beta = beta, log_kappa = log_kappa, log_tau = log_tau,
                         strparam =  strparam,  fields = fields, tmb_silent = tmb_silent,
                         nlminb_silent = nlminb_silent, ...)
    return(res)
}

