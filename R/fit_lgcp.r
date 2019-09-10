#' creates a virtual class that is a superclass to the component classes so then both
#' children inherit from that class
setClassUnion("missing_or_spatialpolygon", c("missing", "SpatialPolygons")) 
#' Function to fit a spatiotemporal log-Gaussian Cox process using TMB and the
#' R_inla namespace for the spde construction of the latent field
#' @param locs a matrix of locations 2xn
#' @param temp.idx numeric aggregated temporal index 1,...,n 
#' @param mesh a Delaney triangulation of the study area created using \code{INLA::inla.mesh.2d()}
#' @param covs a 
#' @param parameters a named vector of parameters. Must include "beta" the regression coefficients of fixed effects,
#' "log_kappa" the smoothness parameter of the random field, if an AR(1) process for time is being fitted then
#' "rho" must also be included as the temporal dependence parameter
#' @param sp  optional spatial polygon of the domain. Should be supplied for lgcp as spde
#' should set weigths to 0 outside of domain
#' @param tmb logical if TRUE will use \code{TMB} framework to fit model (default), otherwise \code{INLA}
#' @param ... arguments to pass into \code{nlminb()} if using TMB, if using INLA framework arguments passed into
#' \code{inla()}
#' @export
setGeneric("fit.lgcp",
           function(locs,temp.idx, mesh, parameters, sp, tmb = TRUE,...){
               standardGeneric("fit.lgcp")
           })

setMethod("fit.lgcp",
          c(locs = "matrix",temp.idx = "numeric",mesh = "inla.mesh",parameters = "vector",
            sp = "missing_or_spatialpolygon" ,
            tmb = "logical"),
          function(locs, temp.idx, mesh, parameters, sp, tmb, ...){
              if(tmb == TRUE){
                  ## check for templates
                  if(!"lgcpar1"%in%getLoadedDLLs()){
                      dll.stelfi()
                  }
                  resp <- list()
                  w.loc <- split(mesh$idx$loc,temp.idx)
                  w.c <- lapply(w.loc,table)
                  for(i in 1:length(w.c)){
                      resp[[i]] <- numeric(mesh$n)
                      count <- as.vector(w.c[[i]])
                      resp[[i]][unique(w.loc[[i]])] <- count
                  }
                  data <- list(resp = resp, ID = as.factor(temp.idx),covariates = covs)
                  spde <- inla.spde2.matern(mesh = mesh,alpha = 2)
                  data$spde <- spde$param.inla[c("M0","M1","M2")]
                  if(class(sp) == "SpatialPolygons"){
                      w <- outwith(mesh = mesh, boundary = sp)
                  }else{
                      w <- rep(1,length(c(diag(data$spde$M0))))
                  }
                  data$area <- w*c(diag(data$spde$M0))
                  ## PARAMS
                  params <- list(beta = parameters[["beta"]],log_kappa = parameters[["log_kappa"]],
                                 x = matrix(0,nrow = mesh$n, ncol = 12),
                                 rho = parameters[["rho"]])
                  fit <- TMB::MakeADFun(data,params,DLL = "lgcpar1",random = c("x"))
                  opt <- stats::nlminb(fit$par,fit$fn,fit$gr, ...)
              }
              return(fit)
          })
              
