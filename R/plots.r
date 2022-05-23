setClassUnion("missing_or_numeric", c("numeric", "missing"))
#' Plots the Hawkes intensty function with decay historical dependence
#' @docType methods
#' @rdname show_hawkes
#' @inheritParams sim_hawkes
#' @param times a vector of numeric observation time points
#' @export
setGeneric("show_hawkes",
           function(times, mu, alpha, beta, marks){
               standardGeneric("show_hawkes")
           })

setMethod("show_hawkes",
          c(times = "vector", mu = "numeric", alpha = "numeric", beta  = "numeric", marks = "vector"),
          function(times, mu, alpha, beta, marks){
              n = length(times)
              max = max(times)
              p = seq(0,max,length.out = 500)
              lam.p = hawke_intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p, marks= marks)
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
show_field <- function(x,smesh, dims = c(500,500), border, colour_option="D",
                       title="Plot"){
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
         ggplot2::scale_fill_viridis_c(option=colour_option) +
         ggplot2::coord_equal()+
         ggplot2::ggtitle(title)
     
     if (!missing(border)){
     plt<- plt+ggplot2::geom_polygon(data=border_f,ggplot2::aes(x=long,y=lat,group=group),fill=NA,color='red')
     }
     plt

 }
#' Plot the lambda estimate from a spatial model

plot_lambda <-function(fit,covariates, smesh,dims=c(500,500), border, colour_option="D",
                       title="Plot"){
        if(!missing(covariates)) {
                designmat <- cbind(1, covariates)
        } else {
                designmat <- matrix(rep(1, smesh$n), ncol = 1)
        }
        
        res = TMB::sdreport(fit)
        field = res$par.random
        beta = res$value[1:ncol(designmat)]
        beta = as.matrix(beta,ncol=1)
        lambda = field + designmat%*%beta
        if (missing(border)){
                stelfi::show_field(x=lambda,smesh=smesh,dims=dims,
                                  colour_option=colour_option, title=title)
        } else {
                stelfi::show_field(lambda,smesh,dims,border,colour_option,title)
        }
        
}
        
        