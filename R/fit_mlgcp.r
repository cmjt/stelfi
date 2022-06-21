#' Fit a marked spatial log-Gaussian Cox process m(LGCP)
#'
#' \code{fit_mlgcp_tmb} fits a marked LGCP using \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field.
#' For a simpler wrapper use \code{\link{fit_mlgcp}} as this is an internal function.
#' @param ypp  The vector of observations.
#' @param marks  A matrix of marks for each observation of the point pattern
#' @param lmat A sparse matrix mapping mesh points to the observations
#' @param strfixed A matrix of fixed structural parameters, defined for each event and mark.
#' Defaults to 1s.
#' Normal distribution: this is the log of standard deviation.
#' Poisson distribution: not used
#' Binomial distribution, this is the number of trials.
#' Gamma distribution, this is the log of the scale.
#' @param methods An integer value:
#' \itemize{
#' \item \code{0} (default), Gaussian distribution, parameter estimated is mean;
#' \item \code{1}, Poisson distribution, parameter estimated is intensity;
#' \item \code{2}, binomial distribution, parameter estimated is logit/probability;
#' \item \code{3}, gamma distribution, the implementation in TMB is shape-scale.
#' }
#' @param betamarks Numeric starting values of the intercept term and
#' covariate coefficients of each mark.
#' @param betapp Numeric starting value of the  intercept of point process.
#' @param marks_coefs_pp The numeric starting value of the coefficient that connects
#' the mark intensity to the point process intensity.
#' @param log_kappa  The numeric starting value of log of kappas for the
#' random field(s). The first element is for the random field of the
#' point process.
#' @param log_tau  The numeric starting value of log of taus for the
#' random field(s). The first element is for the random field of the
#' point process.
#' @param cov_overlap Logical, if \code{TRUE} then  marks and point process
#' share covariates. In that case, To avoid parameter redundancy and non-identifiability
#' {marks_coefs_pp} multiples only the GMRF of the point process.
#' If \code{FALSE}, \code{marks_coef_pp} multiples the full value of lambda
#' @param designmatpp Design matrix for point process. The first column is ones.
#' If there are covariates, then the covariates are in the subsequent columns.
#' @param designmatmarks The design matrix for the marks.
#' @param fields A binary vector indicating whether there is a new random
#' field for each mark. By default, each mark has its own random field.
#' @inheritParams fit_lgcp_tmb
fit_mlgcp_tmb <- function(ypp, marks, lmat, spde, w, strfixed, methods,
                          betamarks, betapp, marks_coefs_pp, log_kappa, log_tau,
                          cov_overlap, designmatpp, designmatmarks, fields,
                          tmb_silent, nlminb_silent, ...) {
    if (!"marked_lgcp" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi("marked_lgcp")
    }
    data <- list(ymarks = marks, ypp = ypp, lmat = lmat,
                 spde = spde$param.inla[c("M0", "M1", "M2")], w = w,
                 methods = methods, designmatpp = designmatpp,
                 designmatmarks = designmatmarks, cov_overlap = cov_overlap,
                 strfixed = strfixed, mark_field = fields)
    param <- list(betamarks = betamarks, betapp = betapp, log_kappa = log_kappa,
                  log_tau = log_tau, marks_coefs_pp = marks_coefs_pp,
                  x = matrix(0, nrow = dim(lmat)[2], ncol = sum(fields) + 1))
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE,
                           random = c("x"), DLL = "marked_lgcp",
                           silent = tmb_silent)
    trace <- if(nlminb_silent) 0 else 1
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                         control = list(trace = trace), ...)
    return(obj)
}
#' Function to fit a spatial  marked
#' log-Gaussian Cox process using \code{TMB}
#'
#' A simple to use wrapper for \code{fit_mlgcp_tmb}
#' @param locs A \code{data.frame} of \code{x} and \code{y} locations, 2xn.
#' @inheritParams fit_lgcp
#' @inheritParams fit_mlgcp_tmb
#' @param parameters a list of named parameters:
#' log_tau, log_kappa, betamarks, betapp, marks_coefs_pp.
#' See \code{\link{fit_mlgcp_tmb}} for further details.
#' @param covariates Covariate(s) corresponding to each area in the spatial mesh
#' @param pp_covariates Which columns of the covariates apply to the point process
#' @param marks_covariates Which columns of the covariates apply to the marks.
#' By default, all covariates apply to the marks only.
#' @inheritParams fit_lgcp_tmb
#' @examples{
#' \dontrun{
#' data(marked, package = "stelfi")
#' loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
#' domain <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(loc.d)),'0')))
#' smesh <- INLA::inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1),
#' max.edge = c(0.3, 0.7), cutoff = 0.05)
#' locs <- cbind(x = marked$x, y = marked$y)
#' marks <- cbind(m1 = marked$m1) ## Gaussian
#' parameters <- list(betamarks = matrix(0, nrow = 1, ncol = ncol(marks)),
#' log_tau = rep(log(1), 2), log_kappa = rep(log(1), 2),
#' marks_coefs_pp = rep(0, ncol(marks)), betapp = 0)
#' fit <- fit_mlgcp(locs = locs, marks = marks,
#' sp = domain, smesh = smesh,
#' parameters = parameters, methods = 0,fields = 1)
#' }
#' }
#' @export
fit_mlgcp <-  function(locs, sp, marks, smesh, parameters=list(), methods,
                       strfixed = matrix(1, nrow = nrow(locs), ncol = ncol(marks)),
                       fields = rep(1, ncol(marks)),
                       covariates, pp_covariates, marks_covariates,
                       tmb_silent = TRUE, nlminb_silent = TRUE, ...) {
    ## Verify args are correct size and class
    n_marks <- ncol(marks)
    n_fields <- sum(fields) + 1
    if (length(log_tau) != n_fields)
        stop("There must be one log_tau for each field")
    if (length(log_kappa) != n_fields)
        stop("There must be one log_kappa for each field")
    if (ncol(betamarks) != n_marks)
         stop("ncol.betamarks must equal ncol.marks")
    if (length(marks_coefs_pp) != n_marks)
        stop("marks_coefs_pp must have length ncol.marks")
    if (length(methods) != n_marks)
        stop("arg methods must have length ncol.marks")
    if (length(fields) != n_marks)
        stop("arg fields must have length ncol.marks")
    if (nrow(strfixed) != nrow(locs))
        stop("nrow.strfixed must be equal to number of points")
    if (ncol(strfixed) != n_marks)
        stop("ncol.strfixed must be equal to ncol.marks")
    if(!missing(covariates)) {
        if(!"matrix" %in% class(covariates))
            stop("arg covariates must be a matrix")
        if (missing(pp_covariates)) {
            pp_covariates <- numeric(length = 0)
        }
        if (missing(marks_covariates)) {
            marks_covariates <- c(1:n_marks)
        }
        if (length(pp_covariates) > ncol(covariates))
            stop("pp_covariates has too many entries")
        if (length(marks_covariates) > ncol(covariates))
            stop("marks_covariates has too many entries")
        if (length(betapp) != (length(pp_covariates) + 1))
            stop("The length of betapp must be one more than the length of pp_covariates")
        if (nrow(betamarks) != (length(marks_covariates) + 1))
            stop("nrow.betamarks must be one more than the length of marks_covariates")
        if(nrow(covariates) != nrow(smesh$loc))
            stop("nrow.covariates must be equal to spatial mesh size")
    } else {
        if(length(betapp) != 1)
            stop("arg betapp must be length 1 if covariates missing")
        if (nrow(betamarks) != 1)
            stop("nrow.betamarks must be 1 if covariates missing")
    }
    
    ## read in parameters
    log_tau <- parameters[["log_tau"]]
    if (is.null(log_tau)) {
      log_tau <- numeric(n_fields)
    }
    log_kappa <- parameters[["log_kappa"]]
    if (is.null(log_kappa)) {
      log_kappa <- numeric(n_fields)
    }
    betamarks <- parameters[["betamarks"]]
    if (is.null(betamarks)) {
      if (!missing(covariates)) {
        betamarks <- matrix(0, nrow = (length(marks_covariates) + 1), ncol = n_marks)
      } else {
        betamarks <- matrix(0, nrow = 1, ncol = n_marks)
      }
    }
    betapp <- parameters[["betapp"]]
    if (is.null(betapp)) {
      area <- sum(get_weights(smesh,sp)$w)
      avg_rate <- log(nrow(locs)/area)
      if (!missing(covariates)) {
        betapp <- numeric(length(pp_covariates))
        betapp[1] <- avg_rate
      } else {
        betapp <- avg_rate
      }
    }
    marks_coefs_pp <- parameters[["marks_coefs_pp"]]
    if (is.null(marks_coefs_pp)) {
      marks_coefs_pp = numeric(n_marks)
    }

    ## data
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys

    ypp <- points.in.mesh(locs, polys)
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    lmat <- INLA::inla.spde.make.A(smesh, locs)
    if(!missing(covariates)) {
        ## overlap of covariates
        if (length(Reduce(intersect, list(marks_covariates, pp_covariates))) > 0) {
          cov_overlap <- 1
        } else {
        cov_overlap <- 0
        }    
        ## Design matrices
        designmatpp <- cbind(1, covariates[,  pp_covariates])
        designmatmarks <- cbind(1, covariates[, marks_covariates])
    } else {
        cov_overlap <- 0
        designmatpp <- matrix(rep(1, length(ypp)), ncol = 1)
        designmatmarks <- matrix(rep(1, length(ypp)), ncol = 1)
    }
    ## Model fitting
    res <- fit_mlgcp_tmb(ypp = ypp, marks = marks, lmat = lmat,
                         spde = spde, w = w_areas, strfixed = strfixed,
                         methods = methods, betamarks = betamarks,
                         betapp = betapp, marks_coefs_pp = marks_coefs_pp,
                         log_kappa = log_kappa, log_tau = log_tau,
                         designmatpp = designmatpp,
                         designmatmarks = designmatmarks, cov_overlap = cov_overlap,
                         fields = fields, tmb_silent = tmb_silent,
                         nlminb_silent = nlminb_silent, ...)
    return(res)
}
