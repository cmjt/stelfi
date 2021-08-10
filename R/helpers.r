setClassUnion("numeric_or_missing", c("numeric", "missing"))
#' Hawkes intensty function with decay historical dependence
#' @inheritParams sim_hawkes
#' @inheritParams show_hawkes
#' @param p vector of pseudo times at which to calculate the intensity
#' @export
setGeneric("hawke_intensity",
           function(mu, alpha, beta, times, p){
           })

setMethod("hawke_intensity",
          c(mu = "numeric",alpha = "numeric" ,beta  = "numeric",times = "vector",
            p = "numeric_or_missing"),
          function(mu, alpha, beta, times, p){
              if(missing(p)) p <- times
              lam <- function(p){
                  mu + alpha*sum(exp(-beta*(p - times))[times < p])
              }
              lam.p <- rep(0,length(p))
              for(i in 1:length(p)){
                  lam.p[i] <- lam(p[i])
              }
              return(lam.p)
          })
#' Coefficient extraction
#' @param object result of a call to \code{fit_hawkes}
#' @export
setGeneric("get_coefs",
           function(object){
               standardGeneric("get_coefs")
           }
           )
setMethod("get_coefs",
          signature(object = "list"),
          function(object){
              summary(TMB::sdreport(object),"report")
          }
          )
