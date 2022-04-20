setClassUnion("numeric_or_missing", c("numeric", "missing"))
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim_hawkes
#' @inheritParams show_hawkes
#' @param p vector of pseudo times at which to calculate the intensity
#' @export
setGeneric("hawke_intensity",
           function(mu, alpha, beta, times, p) {
           })

setMethod("hawke_intensity",
          c(mu = "numeric", alpha = "numeric", beta  = "numeric",
            times = "vector",
            p = "numeric_or_missing"),
          function(mu, alpha, beta, times, p) {
              if (missing(p)) p <- times
              lam <- function(p) {
                  mu + alpha * sum(exp(-beta * (p - times))[times < p])
              }
              lam_p <- rep(0, length(p))
              for (i in seq_along(p)) {
                  lam_p[i] <- lam(p[i])
              }
              return(lam_p)
          })
#' Reported parameter estimates
#' @param object result of a call to \code{fit_hawkes()} or \code{fit_lgcp()}
#' @export
setGeneric("get_coefs",
           function(object) {
               standardGeneric("get_coefs")
           }
           )
setMethod("get_coefs",
          signature(object = "list"),
          function(object) {
              summary(TMB::sdreport(object), "report")
          }
          )
#' Estimated random field(s)
#' @param object result of a call \code{fit_lgcp()}
#' @param plot logical, if TRUE then field(s) plotted
#' @param sd logical, if TRUE then standard errors of field returned
#' @inheritParams fit_lgcp 
#' @export
get_fields <- function(object, smesh, tmesh, plot = FALSE, sd = FALSE) {
    idx <- ifelse(sd, 2, 1)
    x <- summary(TMB::sdreport(object),"random")[,idx]
    if(!missing(tmesh)){
        ind <- rep(seq(tmesh$n), each = smesh$n)
        x <- split(x,ind)
        if(plot) {
            for(i in seq(tmesh$n)) {
                dev.new()
                print(show_field(x[[i]], smesh))
            }
        }
    }else{
        if(plot) print(show_field(x, smesh))
    }
    return(x)
}
      
#' Function to find areas (weights) around the mesh nodes which are
#' within the specified spatial polygon.
#' Relies on inla.mesh.dual function from INLA spde-tutorial
#' \href{SPDE gitbook}{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @export
get_weights <- function(mesh, sp, plot = FALSE){
    dmesh <- inla.mesh.dual(mesh)
    sp::proj4string(dmesh) <- sp::proj4string(sp)
    w <-  sapply(1:length(dmesh), function(i) {
        if (rgeos::gIntersects(dmesh[i, ], sp))
            return(rgeos::gArea(rgeos::gIntersection(dmesh[i, ], sp)))
        else return(0)
    })
    if(missing(plot)) plot = FALSE
    if(plot){
        sp::plot(dmesh, col = "grey")
        sp::plot(mesh, add = TRUE, edge.color = "white")
        sp::plot(sp,add = TRUE)
        points(mesh$loc,pch = 18)
        points(mesh$loc[unlist(w) == 0,],col = "white", pch = 18)
    }
    return(list(weights = unlist(w), polys = dmesh))
}
#' Function that takes in a named  matrix of covariates with
#' ncol equal to the number of covariates
#' returns a list containing the effects ready to be
#' read by \code{inla.stack} and a covariate formula,
#' ready to be read by a call to \code{fit_lgcp_inla}
make.covs <- function(covariates){
    n.covs <- ncol(covariates)
    for(i in 1:n.covs){
        assign(colnames(covariates)[i],covariates[,i],envir = .GlobalEnv)
    }
    cov.effects <- sapply(colnames(covariates),get,simplify = FALSE)
    cov.form <- paste(colnames(covariates),collapse = " + ")
    return(list(cov.effects,cov.form))
}

#' Function that takes in a list of points 
#' (and optionally a weight or covariate for each point)
#' and a mesh of polygons
#' Returns either the number of points in each polygon or sum of the weights

points.in.mesh <- function(xy, dmesh, weights){
  if (missing(weights)){
    sapply(1:length(dmesh), function(i){
      coord <- raster::geom(dmesh[i, ])[,c("x", "y")]
      sum(sp::point.in.polygon(xy[, 1], xy[, 2], coord[, 1], coord[, 2]) > 0)
    })
  }
  else {
    result = rep(0, length(dmesh))
    for (i in 1:length(dmesh)) {
      coord <- raster::geom(dmesh[i, ])[,c("x", "y")]
      temp_result = sp::point.in.polygon(xy[, 1], xy[, 2], coord[, 1], coord[, 2]) > 0
      temp_result = temp_result * weights
      result[i] = sum(temp_result)
    }
    return(result)
  }
}
