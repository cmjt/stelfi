setClassUnion("logical_or_missing", c("logical", "missing")) 
#' Plots the Hawkes intensty function with decay historical dependence
#' @docType methods
#' @rdname show_hawkes
#' @inheritParams sim_hawkes
#' @param times a vector of numeric observation time points
#' @export
setGeneric("show_hawkes",
           function(times,mu, alpha, beta){
               standardGeneric("show_hawkes")
           })

setMethod("show_hawkes",
          c(times = "vector", mu = "numeric", alpha = "numeric", beta  = "numeric"),
          function(times, mu, alpha, beta){
              n = length(times)
              max = max(times)
              p = seq(0,max,length.out = 500)
              lam.p = hawke_intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p)
              ylab = expression(lambda(t))
              col = 1
              lmax = max(lam.p)
              lmin = min(lam.p)
              plot(times,rep(lmin-1,n),ylim = c(lmin-2,lmax),xlab="time",ylab = ylab,col = col,pch = 20)
              lines(p,lam.p,col = "grey")
          })
#' Function to plot or create raster of a given field defined over a mesh
#' @docType methods
#' @rdname show_field
#' @param x vector of length \link{mesh$n} of values at each \link{mesh} node
#' @param dims vector of length 2 defining how fine a projection; default c(300,300)
#' @param col colours of plot; default terrain.colors(100)
#' @inheritParams fit_lgcp_inla
#' @param rast logical; if TRUE create raster object; default FALSE
#' @param legend logical; plot legend; default TRUE
#' @param legend.only logical; legend only to be plotted; default FALSE
#' @param ... arguments to pass into raster (if rast = TRUE) or image()
#' @export
#' @importFrom spatstat.geom as.owin inside.owin
setGeneric("show_field",
           function(x, mesh, dims, col,
                    sp , rast, legend, legend.only, ...){
               standardGeneric("show_field")
           })
setMethod("show_field",
          c(x = "numeric", mesh = "inla.mesh", dims = "numeric", col = "character",
            sp = "missing_or_spatialpolygon",rast = "logical_or_missing",
            legend = "logical_or_missing",legend.only = "logical_or_missing"),
          function(x, mesh, dims, col, sp ,rast, legend, legend.only,...){
              if(missing(rast)) rast = FALSE
              if(missing(legend)) legend = TRUE
              if(missing(legend.only)) legend.only = FALSE
              stopifnot(length(x) == mesh$n)
              proj = INLA::inla.mesh.projector(mesh, dims = dims)
              field.proj = INLA::inla.mesh.project(proj, x)
              if(!missing(sp)){
                  require(maptools)
                  e = expand.grid(proj$x,proj$y)
                  ins =  inside.owin(e[,1],e[,2],as.owin(sp))
                  ins = matrix(ins,nrow=length(proj$x))
                  field.proj[!ins] = NA
              }                                
              if(rast){
                  raster::raster(list(x = proj$x, y=proj$y, z = field.proj),...)
              }else{
                  if(legend.only){
                      fields::image.plot(list(x = proj$x, y=proj$y, z = field.proj),
                                         legend.only = TRUE, col = col,add = TRUE,legend.width = 4,
                                         legend.mar = 0)
                  }else{
                      image(list(x = proj$x, y=proj$y, z = field.proj),col = col, axes = FALSE, ...)
                      if(legend){fields::image.plot(list(x = proj$x, y=proj$y, z = field.proj),
                                                    legend.only = TRUE, col = col,
                                                    legend.shrink = 0.5)}
                      if(!missing(sp)){sp::plot(sp, add = TRUE)}
                  }
              }
          })
