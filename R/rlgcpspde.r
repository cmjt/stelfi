setClassUnion("missingOrnumeric", c("missing", "numeric"))
setClassUnion("SpatialPolygonsDataFrameOrowin", c("SpatialPolygonsDataFrame", "owin")) 
#' Function to simulate from a spatial/spatio temporal 2D LGCP (SPDE model)
#'
#' @return A named matrix (or a list of matricies if spatio-temporal) of point locations and
#' (if a marked point pattern is simulated) a mark values
#' @param sp the \code{SpatialPolygonsDataFrame} or \code{spatstat::owin}
#' for the domain is used to construct the delauney traingulation
#' @param mu numeric, the intercept term to simulate a LGCP, by default is 3
#' @param kappa a numeric constant, parameter of the SPDE model, by deafault is 1
#' @param sigma2 a numeric constant, parameter of the SPDE model, by default this is 1
#' @param n a numeric constant defining the number of time points, by default 1
#' @param rho the AR(1) correlation coefficient for spatio-temporal samples,
#' by default this is 0.9.
#' @export
#' @importFrom spatstat as.owin rLGCP rpoispp
setGeneric("rlgcpspde",
           function(sp,  mu, kappa, sigma2, n, nu, rho){
               standardGeneric("rlgcpspde")
           })
setMethod("rlgcpspde",
          c(sp = "SpatialPolygonsDataFrameOrowin", 
            mu = "missingOrnumeric", kappa = "missingOrnumeric",
            sigma2 = "missingOrnumeric",n = "missingOrnumeric",
            nu = "missingOrnumeric", rho = "missingOrnumeric"),
          function(sp,  mu, kappa, sigma2, n, nu, rho){
              if(missing(mu)) mu <- 3
              if(missing(kappa)) kappa <- 1
              if(missing(sigma2)) sigma2 <- 1
              if(missing(nu)) nu <- 1
              if(missing(n)) n <- 1
              if(class(sp) == "SpatialPolygonsDataFrameOrowin") sp <- as.owin(sp)
              if(n == 1){
                  lg.s <- rLGCP('matern', mu, var = sigma2,
                            scale = kappa / sqrt(8), nu = nu, win = sp)
                  locs <- cbind(x = lg.s$x, y = lg.s$y)
                  return(locs)
              }else{
                  if(missing(rho)) rho <- 0.9
                  lg.s <- rLGCP('matern', mu, var = sigma2, 
                            scale = kappa / sqrt(8), nu = nu, win = sp, nsim = n)
                  vields <- lapply(lg.s,function(x) attributes(x)$Lambda)
                  x <- vields
                  for (j in 2:n){
                      x[[j]] <- rho*x[[j-1]] + sqrt(1-rho^2)*vields[[j]]
                  }
                  pps <- lapply(x,rpoispp)
                  locs <- cbind(x = unlist(sapply(pps, function(x) x$x)),
                                y = unlist(sapply(pps, function(x) x$y)),
                                t = rep(1:n, times = sapply(pps, function(x) x$n)))
                  rownames(locs) <- NULL
                  return(locs)
              }
          }
          )
