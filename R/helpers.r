setClassUnion("numeric_or_missing", c("numeric", "missing"))
#' Hawkes intensity function
#' @inheritParams sim_hawkes
#' @inheritParams show_hawkes
#' @param p An optional vector of pseudo times at which to calculate the intensity.
#' @export
setGeneric("hawke_intensity",
           function(times, mu, alpha, beta, p, marks, background_param) {
           })

setMethod("hawke_intensity",
          c(times = "vector", alpha = "numeric",
            beta  = "numeric", p = "numeric_or_missing"),
          function(times, mu, alpha, beta, 
                   p, marks, background_param) {
              if (missing(p)) p <- times
              lam <- function(p,mark,mu) {
                  mu + mark * alpha * sum(exp(-beta * (p - times))[times < p])
              }
              lam_p <- rep(0, length(p))
              if (class(mu) == "function"){
                mus <- mu(background_param, p)
              } else {
                mus <- rep(mu,length(p))
              }
              for (i in seq_along(p)) {
                  lam_p[i] <- lam(p[i],marks[i],mus[i])
              }
              return(lam_p)
          })
#' Reported parameter estimates
#' 
#' \code{get_coefs} returns the parameter estimates for the fitted model.
#' @param object result of a call to \code{\link{fit_hawkes}} or \code{\link{fit_lgcp}}
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
#' 
#' \code{get_fields} returns 
#' @param object The result of a call to \code{\link{fit_lgcp}}.
#' @param plot Logical, if \code{TRUE} then the returned values are plotted.
#' @param sd Logical, if \code{TRUE} then standard errors of field are returned. 
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
#' Relies on the internal \code{\link{inla.mesh.dual}} function.
#' @seealso \url{https://becarioprecario.bitbucket.io/spde-gitbook/}.
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

#' Internal function that takes in a list of points 
#' (and optionally a weight or covariate for each point)
#' and a mesh of polygons.
#' @param xy A data frame of locations.
#' @param dmesh An object returned by \code{\link{inla.mesh.dual}}.
#' @param weights A vector of weights at each node of \code{dmesh}.
#' @returns Either the number of points in each polygon or sum of the weights.
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
