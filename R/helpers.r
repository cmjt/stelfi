#' creates a virtual class that is a superclass to the component classes so then both children inherit from that class
setClassUnion("numeric_or_NULL", c("numeric", "NULL")) 
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim.hawkes
#' @inheritParams show_hawkes
setGeneric("hawke.intensity",
           function(mu, alpha, beta,times,p = NULL){
           })

setMethod("hawke.intensity",
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
setGeneric("get_weights",
           function(mesh, sp, plot = FALSE){
           })

setMethod("get_weights",
          c(mesh = "inla.mesh", sp = "SpatialPolygonsDataFrame", plot = "logical"),
          function(mesh, sp, plot = FALSE){
              dmesh <- inla.mesh.dual(mesh)
              proj4string(dmesh) <- proj4string(sp)
              w <- sapply(1:length(dmesh), function(i) {
                  if (rgeos::gIntersects(dmesh[i,], sp)){
                      return(as.numeric(sum(sf::st_area(sf::st_intersection(sf::st_as_sf(dmesh[i, ]),
                                                                            lwgeom::st_make_valid(sf::st_as_sf(sp)))))))
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
