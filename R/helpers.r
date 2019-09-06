#' creates a virtual class that is a superclass to the component classes so then both children inherit from that class
setClassUnion("numericOrNULL", c("numeric", "NULL")) 
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim.hawkes
#' @inheritParams plot_hawkes
setGeneric("hawke.intensity",
           function(mu, alpha, beta,times,p = NULL){
           })

setMethod("hawke.intensity",
          c(mu = "numeric",alpha = "numeric" ,beta  = "numeric",times = "vector",p = "numericOrNULL"),
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
