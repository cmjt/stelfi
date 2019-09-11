#' creates a virtual class that is a superclass to the component classes so then both children inherit from that class
setClassUnion("numeric_or_NULL", c("numeric", "NULL")) 
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim.hawkes
#' @inheritParams plot_hawkes
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
#' Function copied from http://inla.r-inla-download.org/r-inla.org/tutorials/spde/R/spde-tutorial-functions.R
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
#' Function to find mesh nodes outwith some spatial polygon
#' relies on inla.mesh.dual function from INLA spde-tutorial
outwith <- function(mesh = NULL,boundary = NULL){
    dmesh <- inla.mesh.dual(mesh)
    proj4string(dmesh) <- proj4string(boundary)
    w <- sapply(1:length(dmesh), function(i) {
        if (rgeos::gIntersects(dmesh[i,], boundary))
            return(rgeos::gArea(rgeos::gIntersection(dmesh[i,], boundary)))
        else return(0)
    })
    return(w)
}
