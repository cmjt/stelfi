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
    x <- summary(TMB::sdreport(fit),"random")[,idx]
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
      
