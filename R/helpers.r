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
            sp::Polygon(p[order(atan2(yy,xx)), ])
        })
        return(sp::SpatialPolygons(lapply(1:mesh$n, function(i)
            sp::Polygons(list(pls[[i]]), i))))
    }
    else stop("It only works for R2!")
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
