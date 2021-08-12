## Needed for registerin S4 methods
setClass("inla.mesh.1d")
setClass("inla.mesh")
setClass("inla.spde")
setClassUnion("spdf_or_sp", c("SpatialPolygonsDataFrame", "SpatialPolygons"))
setClassUnion("numeric_or_missing", c("numeric", "missing"))
setClassUnion("numeric_or_vector", c("numeric", "vector"))
setClassUnion("vector_or_matrix", c("vector", "matrix"))
setClassUnion("inla.mesh.1d_or_missing", c("inla.mesh.1d", "missing"))



#' Function to fit a spatial or spatiotemporal log-Gaussian Cox process using
#' \code{TMB} and the
#' \code{R_inla} namespace for the spde construction of the latent field. For
#' a simpler wrapprt use \code{fit_lgcp()}
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
#' (same length as \code{nrow(designmat)}
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
#' @export
setGeneric("fit_lgcp_tmb",
           function(y, A, designmat, spde, w, idx, beta, x, log_tau,
                    log_kappa, atanh_rho, ...) {
               standardGeneric("fit_lgcp_tmb")
           })
setMethod("fit_lgcp_tmb",
          c(y = "vector_or_matrix", A = "matrix", designmat = "matrix",
            spde = "inla.spde", w = "vector", idx = "vector", beta = "vector",
            x = "matrix", log_tau = "numeric", log_kappa = "numeric",
            atanh_rho = "numeric_or_missing"),
          function(y, A, designmat, spde, w, idx, beta, x, log_tau, log_kappa,
                   atanh_rho, ...) {
              if (!"lgcp" %in% getLoadedDLLs()) {
                  dll_stelfi()
              }
              data <- list(y = y, A = A, designmat = designmat,
                          spde = spde, w = w,
                          idx = idx)
              param <- list(beta = beta, log_tau = log_tau,
                            log_kappa = log_kappa, x = x)
              obj <- TMB::MakeADFun(data = data, parameters = param,
                                    random = c("x"))
              obj$hessian <- TRUE
              opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
              return(obj)
          })
#' Function to fit a spatial or spatiotemporal log-Gaussian Cox process using \code{TMB}
#'
#' A simple to use wrapper for \code{fit_lgcp_tmb}}
#' @param sp \code{SpatialPolygons} of the domain
#' @param locs 2xn \code{data.frame} of locations x, y. If locations have
#' time stapms then this is the third column of the 3xn matrix
#' @param smesh spatial mesh
#' @param tmesh optional, temporal mesh
#' @inheritParams fit_lgcp_tmb
#' @export
setGeneric("fit_lgcp",
           function(locs, sp, smesh, tmesh, beta, log_tau,
                    log_kappa, atanh_rho, ...) {
               standardGeneric("fit_lgcp")
           })
setMethod("fit_lgcp",
          c(locs = "data.frame", sp  = "spdf_or_sp",
            smesh = "inla.mesh",  tmesh = "inla.mesh.1d_or_missing",
            beta = "numeric_or_vector", log_tau = "numeric",
            log_kappa = "numeric",
            atanh_rho = "numeric_or_missing"),
          function(locs, sp, smesh, tmesh, beta, log_tau, log_kappa, atanh_rho, ...) {
               if (!"lgcp" %in% getLoadedDLLs()) {
                  dll_stelfi()
              }
              if (!missing(tmesh)) {
                  stk <- prep(locs = locs, sp = sp, smesh = smesh, tmesh = tmesh)
                   k <- length(tmesh$loc)
              }else{
                  stk <- prep(locs = locs, sp = sp, smesh = smesh)
                  k <- 1
              }
              ## SPDE
              spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
              data <- list(y = stk$data$data$y, A = stk$A,
                                  w = stk$data$data$exposure,
                                  idx = rep(1, length(stk$data$data$y)),
                                  designmat = matrix(1, nrow = length(stk$data$data$y), ncol = 1),
                           spde = spde$param.inla[c("M0", "M1", "M2")])
              param <- list(beta = beta, log_tau = log_tau,
                            log_kappa = log_kappa,  atanh_rho = atanh_rho,
                            x = matrix(0, nrow = spde$n.spde, ncol = k))
              obj <- TMB::MakeADFun(data = data, parameters = param,
                                    random = c("x"), DLL = "lgcp")
              obj$hessian <- TRUE
              opt <- stats::nlminb(obj$par, obj$fn, obj$gr, ...)
              return(obj)

          })

#' Function to prep INLA stuff
#' @inheritParams fit_lgcp
prep <- function(locs, sp, smesh, tmesh) {
    ## E
    dd <- deldir::deldir(smesh$loc[, 1], smesh$loc[, 2])
    tiles <- deldir::tile.list(dd)
    polys <- sp::SpatialPolygons(lapply(seq(1, length(tiles)), function(i) {
        p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
        n <- nrow(p)
        sp::Polygons(list(sp::Polygon(p[c(1:n, 1), ])), i)
    }))
    area <- factor(over(sp::SpatialPoints(cbind(locs$x, locs$y)), polys),
                   levels = seq(1, length(polys)))
    w_areas <- sapply(seq(1, length(tiles)), function(i) {
        p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
        n <- nrow(p)
        pl <- sp::SpatialPolygons(
                      list(sp::Polygons(list(sp::Polygon(p[c(1:n, 1), ])), i)))
        if (rgeos::gIntersects(pl, sp))
            return(rgeos::gArea(rgeos::gIntersection(pl, sp)))
        else return(0)
    })
    ## SPDE
    spde <- INLA::inla.spde2.matern(smesh, alpha = 2)
    ## spatial or spatiotemporal
    if (!missing(tmesh)) {
        k <- length(tmesh$loc)
        idx <- INLA::inla.spde.make.index("s", spde$n.spde, n.group = k)
        w_t <- diag(INLA::inla.mesh.fem(tmesh)$c0)
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
    stk <- inla.stack(
        data = list(y = agg_dat$Freq, exposure = e0),
        A = list(1, a_st),
        effects = list(idx, list(b0 = rep(1, nrow(agg_dat)))))
    return(stk)
}
