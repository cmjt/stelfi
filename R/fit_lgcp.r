#' Creates a virtual class that is a superclass to the component classes so then both
#' children inherit from that class
setClassUnion("missing_or_spatialpolygon", c("missing", "SpatialPolygonsDataFrame","SpatialPolygons"))
setClassUnion("null_or_factor", c("NULL", "factor"))
# Needed for registerin S4 methods
setClass("inla.mesh")
#' Function to fit a spatiotemporal log-Gaussian Cox process using TMB and the
#' R_inla namespace for the spde construction of the latent field
#' @param temp.idx numeric aggregated temporal index 1,...,n
#' @inheritParams fit_lgcp_inla
#' @param covs a list of covariates at the mesh nodes for each time index
#' @param parameters a named list of parameters. Must include \code{beta} the regression coefficients of fixed effects,
#' "log_kappa" the smoothness parameter of the random field, if an AR(1) process for time is being fitted then
#' "rho" must also be included as the temporal dependence parameter
#' @param sp  optional spatial polygon of the domain
#' should set weigths to 0 outside of domain
#' @param ... arguments to pass into \code{nlminb()}
#' @export
setGeneric("fit_lgcp_tmb",
           function(locs,temp.idx, mesh, parameters, covs, sp, ...){
               standardGeneric("fit_lgcp_tmb")
           })

setMethod("fit_lgcp_tmb",
          c(locs = "matrix",temp.idx = "null_or_factor",mesh = "inla.mesh",parameters = "list",
            covs = "list",sp = "missing_or_spatialpolygon"),
          function(locs, temp.idx, mesh, parameters, covs, sp, ...){
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
              if(class(sp) == "SpatialPolygonsDataFrame"){
                  w <- get_weights(mesh = mesh,sp,FALSE)
                  data$area <- w*c(Matrix::diag(data$spde$M0))
              }else{
                  data$area <- c(Matrix::diag(data$spde$M0))
              }
              ## params
              params <- list(beta = parameters[["beta"]],log_kappa = parameters[["log_kappa"]],
                             x = matrix(0,nrow = mesh$n, ncol = length(table(temp.idx))),
                             rho = parameters[["rho"]])
              fit <- TMB::MakeADFun(data,params,DLL = "lgcpar1",random = c("x"))
              opt <- stats::nlminb(fit$par,fit$fn,fit$gr, ...)
              
              return(fit)
          })

#' Function to simulate a spatiotemporal log-Gaussian Cox process using TMB and the
#' R_inla namespace for the spde construction of the latent field from a fitted model
#' or from a template
#' @param x a fitted model from a call to \code{fit_lgcp_tmb} 
#' @inheritParams fit_lgcp_tmb
#' @export
setGeneric("sim_lgcp",
           function(x,locs, temp.idx, mesh, parameters, covs, sp){
               standardGeneric("sim_lgcp")
           })
setMethod("sim_lgcp",
          c(x = "list",locs = "missing", temp.idx = "missing", mesh = "missing", parameters = "missing",covs = "missing",
            sp = "missing"),
          function(x, locs, temp.idx, mesh, parameters, covs, sp){
              simdata = x$simulate(par = x$obj,complete = TRUE)
              return(simdata)
          })
setMethod("sim_lgcp",
          c(x = "list", locs = "missing", temp.idx = "missing", mesh = "missing", parameters = "list", covs = "missing",
            sp = "missing"),
          function(x, locs, temp.idx, mesh, parameters, covs, sp){
              simdata = x$simulate(par = parameters,complete = TRUE)
              return(simdata)
          })
setMethod("sim_lgcp",
          c(x = "missing", locs = "matrix",temp.idx = "factor",mesh = "inla.mesh",parameters = "list",covs = "list",
            sp = "missing_or_spatialpolygon"),
          function(x, locs, temp.idx, mesh, parameters, covs, sp){
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
              if(class(sp) == "SpatialPolygonsDataFrame"){
                  w <- outwith(mesh = mesh, boundary = sp)
                  data$area <- w*c(Matrix::diag(data$spde$M0))
              }else{
                  data$area <- c(Matrix::diag(data$spde$M0))
              }
              ## params
              params <- list(beta = parameters[["beta"]],log_kappa = parameters[["log_kappa"]],
                             x = matrix(0,nrow = mesh$n, ncol = length(table(temp.idx))),
                             rho = parameters[["rho"]])
              onj <- TMB::MakeADFun(data,params,DLL = "lgcpar1",random = c("x"))
              simdata = x$simulate(par = parameters,complete = TRUE)
              return(simdata)
          })
#' Wrapper for \code{INLA} to fit a LGCP
#' #' @return A \code{INLA::inla} result object
#'
#' @param mesh an `inla.mesh` object i.e. delauney triangulation of the domain, an
#' object returned by \link{INLA::inla.mesh.2d}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param sp spatial polygon of the point pattern observation window (optional). if supplied
#' weights at the mesh nodes outwith this will be set to zero.
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T) (optional).
#' @param covariates a named data.frame of covariates (optional) 
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{INLA:::pcprior} is used with \code{param = c(0.9,0.9)}.
#' @param prior.range pc prior for the range of the latent field supplied as the vector c(range0,Prange)  (i.e., P(range < range0) = Prange), by default is \code{ c(5,0.9)}. NOTE should be changed to reflect range of the domain. 
#' @param prior.sigma pc prior for the sd of the latent field supplied as a vector (sigma0,Psigma) (i.e., P(sigma > sigma0) = Psigma), by default is \code{ c(1,0.005)}.
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param control.inla a list to control model fitting (as per inla)
#' @param control.fixed a list as per inla by default sets prior for precision intercept
#' @param ... other arguments taken by \code{inla}
#' @importMethodsFrom Matrix diag
#' @export
fit_lgcp_inla <- function(mesh, locs, sp, temp = NULL, covariates = NULL,
                          prior.rho = list(theta = list(prior='pccor1', param = c(0.0, 0.9))),
                          prior.range = c(5,0.9) ,
                          prior.sigma = c(1,0.005),
                          verbose = FALSE,
                          control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                          control.fixed = list(prec.intercept = 0.001),
                          return.attributes = FALSE,
                          ...){
    if(class(sp) == "SpatialPolygons"){
        sp <- SpatialPolygonsDataFrame(sp,data = data.frame(rep(1,length(sp))))
    }
    mesh <- mesh
    spde <- inla.spde2.pcmatern(mesh = mesh,
                                prior.range = prior.range,
                                prior.sigma = prior.sigma)
    ## number of observations
    n <- nrow(locs)
    ## number of mesh nodes
    nv <- mesh$n
    ## weights at mesh nodes
    weights <-  get_weights(mesh, sp, FALSE)
    if(!is.null(temp)){
        k <- max(temp)
        mesh.t <- inla.mesh.1d(seq(1, k, by = 1))
        Ast <- inla.spde.make.A(mesh = mesh,
                                loc = locs,
                                n.group = k,
                                group = temp,
                                group.mesh = mesh.t)
        field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
        st.vol <- rep(weights, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), nv)
        ctr.g <- list(model = 'ar1',param = prior.rho)
        y.pp <- rep(0:1, c(k * nv, n))
        expected <- c(st.vol, rep(0, n))
    }else{
        ## exposure (E)
        expected <- c(weights, rep(0, n))
        ## vector for observations
        y.pp <- rep(0:1, c(nv, n))
        ## integration points
        imat <- Diagonal(nv, rep(1, nv))
        ## projection matrix for observed points
        lmat <- inla.spde.make.A(mesh, locs)
        ## entire projection matrix
        A.pp <- rbind(imat, lmat)
    }
    if(!is.null(covariates)){
        m <- stelfi:::make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y=y.pp, e=expected),
                            A=list(A.pp,1,1),
                            effects=list(field = 1:nv,
                                         b0 = rep(1,length(y.pp)),
                                         cov.effects = cov.effects))
        if(!is.null(temp)){
            formula <- paste("y", "~  0  + b0 +", cov.form ,
                             " + f(field, model = spde, group = field.group,control.group = ctr.g)")
        }else{
            formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
        }
    }else{
        if(!is.null(temp)){
            stack <- inla.stack(
                data = list(y = y.pp, e = expected), 
                A = list(rbind(Diagonal(n = k * nv), Ast), 1), 
                effects = list(field, list(b0 = rep(1, k * nv + n))))
            formula <- y ~ 0 + b0 + f(field, model = spde, group = field.group,
                                      control.group = ctr.g)
        }else{
            stack <- inla.stack(
                data = list(y = y.pp, e = expected), 
                A = list(1, A.pp),
                effects = list(list(b0 = rep(1, nv + nrow(locs))), 
                               list(field = 1:nv)),
                tag = 'pp')
            formula <- y ~ 0  + b0 + f(field, model = spde)
        }
    }
    ##call to inla
    result <- inla(as.formula(formula), family = "poisson",
                   data = inla.stack.data(stack),
                   E = inla.stack.data(stack)$e,
                   control.predictor = list(A = inla.stack.A(stack),compute = TRUE),
                   control.inla = control.inla,
                   control.fixed = control.fixed,
                   verbose = verbose,
                   ...)
    if(return.attributes) {
        attributes(result)$mesh <- as.list(mesh)
        attributes(result)$weights <- weights
    }
    
    result
}
