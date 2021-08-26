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
show_field <- function(x,smesh, dims = c(500,500)){
    nx <- dims[1]
    ny <- dims[2]
    px <- inlabru::pixels(smesh, nx = nx, ny = ny)
    A <- INLA::inla.spde.make.A(smesh, px)
    px$color <- as.vector(A %*% x)
    plt <- ggplot2::ggplot(as.data.frame(px), ggplot2::aes(x, y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = color)) +
        ggplot2::coord_equal()
    plt
}
