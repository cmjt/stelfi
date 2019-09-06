#' Plots the Hawkes intensty function with decay historical dependence
#' @docType methods
#' @rdname plot_hawkes
#' @inheritParams sim.hawkes
#' @param times a vector of numeric observation time points
#' @export
#'
setGeneric("plot_hawkes",
           function(times,mu, alpha, beta){
               standardGeneric("plot_hawkes")
           })

setMethod("plot_hawkes",
          c(times = "vector", mu = "numeric", alpha = "numeric", beta  = "numeric"),
          function(times, mu, alpha, beta){
              n = length(times)
              max = max(times)
              p = seq(0,max,length.out = 500)
              lam.p = hawke.intensity(mu = mu, alpha = alpha, beta = beta, times = times, p = p)
              ylab = expression(lambda(t))
              col = 1
              lmax = max(lam.p)
              lmin = min(lam.p)
              plot(times,rep(lmin-1,n),ylim = c(lmin-2,lmax),xlab="time",ylab = ylab,col = col,pch = 20)
              lines(p,lam.p,col = "grey")
          })
