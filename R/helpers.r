#' creates a virtual class that is a superclass to the component classes so then both children inherit from that class
setClassUnion("numeric_or_NULL", c("numeric", "NULL")) 
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim_hawkes
#' @inheritParams show_hawkes
#' @export
setGeneric("hawke_intensity",
           function(mu, alpha, beta,times,p = NULL){
           })

setMethod("hawke_intensity",
          c(mu = "numeric",alpha = "numeric" ,beta  = "numeric",times = "vector",
            p = "numeric_or_NULL"),
          function(mu, alpha, beta, times, p){
              if(is.null(p)) p <- times
              lam <- function(p){
                  mu + alpha*sum(exp(-beta*(p-times))[times<p])
              }
              lam.p <- rep(0,length(p))
              for(i in 1:length(p)){
                  lam.p[i] <- lam(p[i])
              }
              lam.p
          })
#' Function copied from \url{http://www.r-inla.org/spde-book}
#' which is suppied alongside the
#' \href{SPDE gitbook}{https://becarioprecario.bitbucket.io/spde-gitbook/}
inla.mesh.dual <- function(mesh) {
    if (mesh$manifold=='R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        library(parallel)
        pls <- mclapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[,k]==i)
                if (length(j)>0) 
                    return(rbind(ce[j, , drop=FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                else return(ce[j, , drop=FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[,1]==i)
            j2 <- which(mesh$segm$bnd$idx[,2]==i)
            if ((length(j1)>0) | (length(j2)>0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
            }
            else {
                yy <- p[,2]-mesh$loc[i, 2]
                xx <- p[,1]-mesh$loc[i, 1]
            }
            Polygon(p[order(atan2(yy,xx)), ])
        })
        return(SpatialPolygons(lapply(1:mesh$n, function(i)
            Polygons(list(pls[[i]]), i))))
    }
    else stop("It only works for R2!")
}
#' Function to find areas (weights) around the mesh nodes which are
#' within the specified spatial polygon.
#' Relies on inla.mesh.dual function from INLA spde-tutorial
#' \href{SPDE gitbook}{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @export
#' @importFrom sf st_area st_intersection st_as_sf st_make_valid
#' @importFrom rgeos gIntersects
setGeneric("get_weights",
           function(mesh, sp, plot = FALSE){
           })

setMethod("get_weights",
          c(mesh = "inla.mesh", sp = "SpatialPolygonsDataFrame", plot = "logical"),
          function(mesh, sp, plot = FALSE){
              dmesh <- inla.mesh.dual(mesh)
              proj4string(dmesh) <- proj4string(sp)
              w <- sapply(1:length(dmesh), function(i) {
                  if (gIntersects(dmesh[i,], sp)){
                      return(as.numeric(sum(st_area(st_intersection(st_as_sf(dmesh[i, ]),
                                                                    st_make_valid(st_as_sf(sp)))))))
                  }
                  else return(0)
              })
              if(plot){
                  sp::plot(dmesh, col = "grey")
                  sp::plot(mesh, add = TRUE, edge.color = "white")
                  sp::plot(sp,add = TRUE)
                  points(mesh$loc,pch = 18)
                  points(mesh$loc[unlist(w) == 0,],col = "white", pch = 18)
              }
              return(unlist(w))
          })
#' Function to extract fields from a fitted model (INLA only ATM)
#' @param x an \code{inla} object
#' @param mesh an object of class \code{inla.mesh}
#' @param t optional, if supplied specifies the number of time steps
#' @param mean logical, if TRUE extracts point estimates of the fields,
#' else returns the standard errors
#' @export
setGeneric("get_fields",
           function(x, mesh, t = NULL, mean){
           })
setGeneric("get_fields",
           function(x = "inla",mesh = "inla.mesh", t = "numeric_or_NULL", mean = "logical"){
               field.names <- names(x$summary.random)
               if(mean){
                   fields <- lapply(1:length(field.names),
                                    function(f) x$summary.random[[field.names[f]]]$mean)
               }else{
                   fields <- lapply(1:length(field.names),
                                    function(f) x$summary.random[[field.names[f]]]$sd)
               }
               if(is.numeric(t))  fields <- lapply(fields, split, rep(1:t,each = mesh$n))
               names(fields) <- field.names
               return(fields)
           })
#'  \url{https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html}
owin_to_polygons <- function(x, id = "1") {
    stopifnot(is.owin(x))
    x <- as.polygonal(x)
    closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
    pieces <- lapply(x$bdry,
                     function(p) {
                         Polygon(coords = closering(cbind(p$x,p$y)),
                                 hole = is.hole.xypolygon(p))  })
    z <- Polygons(pieces, id)
    
    return(z)
}
#' Function to convert a\code{spatstat} \code{owin} object to
#' a \code{SpatialPolygonsDataFrame}
#' @source \url{https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html}
#' @param x an object of class \code{owin}
#' @return a \code{SpatialPolygonsDataFrame}
#' @export
owin_to_sp <- function(x) {
    require(spatstat.utils)
    stopifnot(is.owin(x))
    y <- owin_to_polygons(x)
    z <- SpatialPolygonsDataFrame(SpatialPolygons(list(y)),
                                  data = data.frame(rep(1,length(list(y)))))
    return(z)
}
#' Function that takes in a named  matrix of covariates with
#' ncol equal to the number of covariates
#' returns a list containing the effects ready to be
#' read by \code{inla.stack} and a covariate formula,
#' ready to be read by a call to \code{inla}
make.covs <- function(covariates){
    n.covs = ncol(covariates)
    for(i in 1:n.covs){
        assign(colnames(covariates)[i],covariates[,i],envir = .GlobalEnv)
    }
    cov.effects = sapply(colnames(covariates),get,simplify = FALSE)
    cov.form = paste(colnames(covariates),collapse = " + ")
    return(list(cov.effects,cov.form))
}
