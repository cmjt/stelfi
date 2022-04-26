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
#' @param strfixed a matrix of fixed structural parameters.
#' Normal distribution: this is the log of element-wise standard deviation.
#' Poisson distribution: not used
#' Binomial distribution, this is the number of trials.
#' The default input from R should be 1s, i.e., the Bernoulli distribution.
#' Gamma distribution, this is the log of scale.
#' @param methods integer:
#' 0 - normal distribution. Parameter estimated is mean
#' 1 - poisson distribution. Parameter estimated is intensity
#' 2 - binomial distribution. Parameter estimated is logit/probability
#' 3 - gamma distribution, the implementation in TMB is shape-scale. 
#' @param betaresp starting value of the intercept term of each response
#' @param betapp starting value of the  intercept of point process
#' @param beta_coefs_pp slope of mark lambda versus point process lambda
#' @param log_kappa  log of kappas for the random field. The lengths of these vectors are no. of resp + 1.
#' The first element is for the random field of the point process.
#' @param log_tau  log of taus for the random field. The lengths of these vectors are no. of resp + 1.
#' The first element is for the random field of the point process.
#' @param designmat Design matrix. The first column is ones.
#' If there are covariates, then the covariates are in the subsequent columns
#' @param fields Binary vector indicating whether there is a GMRF for each mark. 
#' @param tmb_silent logical, default `TRUE`:
#' TMB inner optimization tracing information will be printed.
#' @param nlminb_silent logical, default `TRUE`:
#' print function and parameters every iteration.
#' @param ... arguments to pass into \code{nlminb()}
#' @export
fit_mlgcp_tmb <- function(ypp, marks, lmat, spde, w, strfixed, methods,
                          betaresp, betapp, beta_coefs_pp, log_kappa, log_tau,
                          designmat, fields, tmb_silent,
                          nlminb_silent, ...){
    if (!"prefsampling" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("prefsampling")
    }
    data <- list(yresp = marks, ypp = ypp, lmat = lmat,
                 spde = spde$param.inla[c("M0", "M1", "M2")], w = w,
                 methods = methods, designmat = designmat,
                 strfixed = strfixed, mark_field = fields)
    param <- list(betaresp = betaresp, betapp = betapp, log_kappa = log_kappa, 
                  log_tau = log_tau, beta_coefs_pp = beta_coefs_pp,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = sum(fields) + 1))
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
#' @param locs 2xn \code{data.frame} of locations x, y. If locations have
#' time stapms then this is the third column of the 3xn matrix
#' @param sp \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} of the domain
#' @param smesh spatial mesh of class \code{"inla.mesh"}
#' @param parameters a list of named parmeters:
#' log_tau, log_kappa, betaresp, betapp, beta_coefs_pp. 
#' See fit_mlgcp_tmb header for further details
#' \code{smesh} and \code{tmesh} node combination.
#' @param fields See fit_mlgcp_tmb header
#' @param covariates Covariate(s) corresponding to each area in the spatial mesh
#' @inheritParams fit_lgcp_tmb
#' @export
fit_mlgcp <-  function(locs, sp, marks, smesh, parameters, methods,
                       strfixed, fields, covariates,
                       tmb_silent = TRUE,
                       nlminb_silent = TRUE, ...) {
    
    ## convert svs
    log_tau <- parameters[["log_tau"]]
    log_kappa <- parameters[["log_kappa"]]
    betaresp <- parameters[["betaresp"]]
    betapp <- parameters[["betapp"]]
    beta_coefs_pp <- parameters[["beta_coefs_pp"]]
    
    ## Verify args are correct size and class
    if(!"matrix" %in% class(idx)) stop("arg idx must be a matrix")
    n_marks <- ncol(marks)
    if (length(log_tau) != (n_marks+1)) stop("log_tau must have length ncol.marks+1")
    if (length(log_kappa) != (n_marks+1)) stop("log_kappa must have length ncol.marks+1")
    if (length(betaresp) != (n_marks)) stop("betaresp must have length ncol.marks")
    if (length(beta_coefs_pp) != (n_marks)) stop("beta_coefs_pp must have length ncol.marks")
    if (length(methods) != n_marks) stop("arg methods must have length ncol.marks")
    if (length(fields) != n_marks) stop("arg fields must have length ncol.marks")
    if (nrow(strfixed) != nrow(locs)) stop("nrow.strfixed must be equal to number of points")
    if (ncol(strfixed) != n_marks) stop("ncol.strfixed must be equal to ncol.marks")
    if(!missing(covariates)){
        if(!"matrix" %in% class(covariates)) stop("arg covariates must be a matrix")
        if(length(betapp) != (ncol(covariates) + 1)) stop("arg betapp should be length ncol.covariates + 1")
        if(nrow(covariates) != nrow(smesh$loc))
            stop("nrow.covariates must be equal to spatial mesh size")
    } else {
        if(length(betapp) != 1) stop("arg betapp should be length 1 if covariates missing")
    }
    
    ## data
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys

    ypp <- stelfi::points.in.mesh(locs, polys)
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    
    ## Designmat
    if(!missing(covariates)) {
        designmat <- cbind(1, covariates)
    } else {
        designmat <- matrix(rep(1, length(ypp) ), ncol = 1)
    }
    
    ## Model fitting
    res <- fit_mlgcp_tmb(ypp = ypp, marks = marks, lmat = lmat,
                         spde = spde,
                         w = w_areas, idx = idx, strfixed = strfixed,
                         methods = methods,
                         betaresp = betaresp, betapp = betapp, beta_coefs_pp = beta_coefs_pp,
                         log_kappa = log_kappa, log_tau = log_tau, designmat = designmat,
                         fields = fields, tmb_silent = tmb_silent,
                         nlminb_silent = nlminb_silent, ...)
    return(res)
}

