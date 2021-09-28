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
#' @param methods integer:  0 - normal distribution, strparam/strfixed as log_sigma.
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
#' @param strparam see \code{strfixed}
#' @param silent logical, by default FALSE. If TRUE model fitting progress
#' not printed to console.
#' @param ... arguments to pass into \code{nlminb()}
#' @export
fit_mlgcp_tmb <- function(ypp, marks, lmat, spde, w, idx, strfixed, methods,
                          betaresp, betapp, beta, log_kappa, log_tau,
                          strparam, silent = FALSE, ...){
    if (!"prefsampling" %in% getLoadedDLLs()) {
        stelfi::dll_stelfi()
    }
    data <- list(yresp = yresps, ypp = ypp, lmat = lmat,
                 spde = spde$param.inla[c("M0", "M1", "M2")], w = w,
                 idx = idx, methods = c(0, 2, 2, 3),
                 strfixed = cbind(rep(log(0.25), dim(yresps)[1]), 1, 1, 2))
    param <- list(betaresp = rep(0, n2), betapp = 0, beta = matrix(0, nrow = n2, ncol = n2 + 1),
                  log_kappa = rep(0, dim(idx)[2]), log_tau = rep(0, dim(idx)[2]), strparam = rep(0, n2),
                  x = matrix(0, nrow = dim(lmat)[2], ncol = sum(diag(idx[, -1]) > 0) + 1))
    obj <- TMB:::MakeADFun(data, param, hessian = TRUE, random = c("x"))
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
#' @param tau tau parameter for the GMRF
#' @param kappa kappa parameter for the GMRF
#' \code{smesh} and \code{tmesh} node combination.
#' @inheritParams fit_lgcp_tmb
#' @export
fit_mlgcp <-  function(locs, sp, marks, smesh, beta, tau, kappa, rho,  silent, ...) {
    if(missing(silent)) silent <- FALSE
    if(sum(names(locs) %in% c("x","y")) < 2) stop("Named variables x and y required in arg locs")
    ## data
    n2 <- dim(marks)[2]
    idx <- cbind(c(1, 1, 0, 1), matrix(0, nrow = n2, ncol = n2))
    idx[3, 4] <- 1
    tmp <- prep_data_marked(locs = locs, sp = sp, smesh = smesh)
    ## SPDE
    stk <- tmp[[1]]
    a_st <- tmp[[2]]
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## convert svs
    log_tau <- log(tau)
    log_kappa <- log(kappa)


    ## Model fitting
    res <- fit_mlgcp_tmb(ypp = stk$data$data$y, marks = marks, lmat = lmat,
                         spde = spde$param.inla[c("M0", "M1", "M2")],
                         w = stk$data$data$exposure, idx = idx, strfixed = strfixed,
                         methods = methods,
                         betaresp = betaresp, betapp = betapp,
                         beta = beta, log_kappa = log_kapp, log_tau = log_tau,
                         strparam = strparam, x = x, silent = silent, ...)
    return(res)
}



#' Function to prep data as per INLA stack
#' @inheritParams fit_mlgcp
prep_data_marked <- function(locs, sp, smesh) {
    ## E
    w <- get_weights(mesh = smesh, sp = sp, plot = FALSE)
    w_areas <- w$weights
    polys <- w$polys
    area <- factor(sp::over(sp::SpatialPoints(cbind(locs$x, locs$y)), polys),
                   levels = seq(1, length(polys)))
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    agg_dat <- as.data.frame(table(area))
    agg_dat[[1]] <- as.integer(as.character(agg_dat[[1]]))
    e0 <- w_areas[agg_dat$area]
    a_st <- INLA::inla.spde.make.A(smesh, smesh$loc[agg_dat$area, ])
    stk <- INLA::inla.stack(
                     data = list(y = agg_dat$Freq, exposure = e0),
                     A = list(1, a_st),
                     effects = list(idx,b0 = rep(1, nrow(agg_dat))))
    return(list(stk,a_st))
}
