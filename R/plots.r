setClassUnion("missing_or_numeric", c("numeric", "missing"))
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
              data = data.frame(x = p, y = lam.p)
              line = ggplot2::ggplot(data = data, ggplot2::aes(x = .data$x, y = .data$y)) +
                  ggplot2::xlab("") +
                  ggplot2::ylab(expression(lambda(t))) + 
                  ggplot2::geom_line() +  ggplot2::theme_minimal()
              hist =  ggplot2::ggplot(data = data.frame(times = times),  ggplot2::aes(x = times)) +
                  ggplot2::geom_histogram() +  ggplot2::theme_minimal() +
                  ggplot2::xlab("times") +  ggplot2::ylab("Number of events")
              gridExtra::grid.arrange(line, hist, ncol = 1)
          })
#' Plots estimated random field(s) of a LGCP
#' @docType methods
#' @rdname show_field
#' @param x a vector of values for each mesh node
#' @param dims by default \code{c(500,500)}, vector of length 2 specifying
#' spatial pixel resolution
#' @inheritParams get_fields
#' @export
show_field <- function(x,smesh, dims = c(500,500), border){
     nx <- dims[1]
     ny <- dims[2]
     px <- inlabru::pixels(smesh, nx = nx, ny = ny)
     A <- INLA::inla.spde.make.A(smesh, px)
     px$color <- as.vector(A %*% x)
     # Convert border data from Spatial Polygons to dataframe
     if (!missing(border)) {
        border_f<-ggplot2::fortify(border)
        }

     plt <- ggplot2::ggplot(as.data.frame(px), ggplot2::aes(x, y)) +
         ggplot2::geom_tile(ggplot2::aes(fill = color)) +
         ggplot2::scale_fill_viridis_c() +
         ggplot2::coord_equal()
     
     if (!missing(border)){
     plt<- plt+ggplot2::geom_polygon(data=border_f,ggplot2::aes(x=long,y=lat,group=group),fill=NA,color='red')
     }
     plt

 }

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
# #setGeneric("show_field",
#            function(x, mesh, dims, col,
#                     sp , rast, legend, legend.only, ...){
#                standardGeneric("show_field")
#            })
# #setMethod("show_field",
#           c(x = "numeric", mesh = "inla.mesh", dims = "numeric", col = "character",
#             sp = "missing_or_spatialpolygon",rast = "logical_or_missing",
#             legend = "logical_or_missing",legend.only = "logical_or_missing"),
#show_field <- function(x, mesh, dims, col, sp ,rast, legend, legend.only,...){
#              if(missing(rast)) rast = FALSE
#              if(missing(legend)) legend = TRUE
#              if(missing(legend.only)) legend.only = FALSE
#              stopifnot(length(x) == mesh$n)
#              proj = INLA::inla.mesh.projector(mesh, dims = dims)
#              field.proj = INLA::inla.mesh.project(proj, x)
#              if(!missing(sp)){
#                  require(maptools)
#                  e = expand.grid(proj$x,proj$y)
#                  ins =  inside.owin(e[,1],e[,2],as.owin(sp))
#                  ins = matrix(ins,nrow=length(proj$x))
#                  field.proj[!ins] = NA
#              }                                
#              if(rast){
#                  raster::raster(list(x = proj$x, y=proj$y, z = field.proj),...)
#              }else{
#                  if(legend.only){
#                      fields::image.plot(list(x = proj$x, y=proj$y, z = field.proj),
#                                         legend.only = TRUE, col = col,add = TRUE,legend.width = 4,
#                                         legend.mar = 0)
#                  }else{
#                      image(list(x = proj$x, y=proj$y, z = field.proj),col = col, axes = FALSE, ...)
#                      if(legend){fields::image.plot(list(x = proj$x, y=proj$y, z = field.proj),
#                                                    legend.only = TRUE, col = col,
#                                                    legend.shrink = 0.5)}
#                      if(!missing(sp)){sp::plot(sp, add = TRUE)}
#                  }
#              }
#          }
