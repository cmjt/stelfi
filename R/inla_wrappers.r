#' Function to fit a  spatio-temporal log-Gaussian Cox process model  with an intercept and covariates (optional)
#'
#' @return A \code{inla} result object
#' @inheritParams fit_lgcp
#' @param locs a named ("x", and "y") data frame of observation locations,
#' where each row corresponds to the observation. If temporal dimension is included
#' "t" sould be the additional column.
#' @param covariates a named data.frame of covariates 
#' @param prior.rho prior for the temporal correlation coefficient,
#' by default a \code{pcprior} is used with \code{param=c(0,0.9)}.
#' @param prior.range pc prior for the range of the latent field (rnage0,Prange)
#' (i.e., P(range < range0) = Prange NOTE should be changed to reflect range
#' of the domain by default this is 400m 
#' @param prior.sigma pc prior for the sd of the latent field (sigma0,Psigma)
#' by default c(1,0.05) i.e., prob sigma > 1 = 0.05 
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param control.inla a list to control model fitting (as per inla)
#' @param control.fixed a list as per inla by default sets prior for precision intercept
#' @param ... other arguments taken by inla
#' @importMethodsFrom Matrix diag
#' @export

fit_lgcp_inla <- function(smesh = NULL, locs = NULL, domain = NULL,
                          covariates = NULL, tmesh = NULL,
                          prior.rho = list(theta = list(prior='pc.cor1', param = c(0.7, 0.7))),
                          prior.range = c(0.05,0.01) ,
                          prior.sigma = c(1,0.01),
                          verbose = FALSE,
                          control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                          control.fixed = list(prec.intercept = 0.001),
                          ...){
    spde <- inla.spde2.pcmatern(mesh = smesh,
                                prior.range = prior.range,
                                prior.sigma = prior.sigma) 
    ## number of observations
    n <- nrow(locs)
    xy <- as.matrix(locs[, c("x", "y")])
    ## number of spatial mesh nodes
    nv <- smesh$n
    if(!is.null(tmesh)){
        k <- tmesh$n
        temp <- locs[, "t"]
        A.pp <- inla.spde.make.A(mesh = smesh, loc = xy, n.group = k,
                                 group = temp, group.mesh = tmesh)
        idx <- inla.spde.make.index('field', n.spde = spde$n.spde, group = temp, n.group = k)
        w <- stelfi::get_weights(smesh, domain)$weights
        st.vol <- rep(w, k) * rep(Matrix::diag(inla.mesh.fem(mesh.t)$c0), nv)
        expected <- c(st.vol, rep(0, n))
        ctr.g <- list(model = 'ar1',param = prior.rho)
        y.pp <- rep(0:1, c(k * nv, n))
    } else {
        w <- stelfi::get_weights(smesh, domain)$weights
        y.pp <- rep(0:1, c(nv, n))
        expected <- c(w, rep(0, n))
        imat <- Diagonal(nv, rep(1, nv))
        lmat <- inla.spde.make.A(mesh, xy)
        A.pp <- rbind(imat, lmat)
    }
    if(!is.null(covariates)){
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data = list(y = y.pp, e = expected),
                            A=list(A.pp,1,1),
                            effects = list(field = field, b0 = rep(1,length(y.pp)),
                                           cov.effects = cov.effects))
        if(!is.null(tmesh)){
            formula <- paste("y", "~  0  + b0 +", cov.form ,
                             " + f(field, model = spde, group = field.group, control.group = ctr.g)")
        }else{
            formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
        }
    }else{
        if(!is.null(tmesh)){
            stack <- inla.stack(
                data = list(y = y.pp, e = expected), 
                A = list(rbind(Diagonal(n = k * nv), A.pp), 1), 
                effects = list(idx, list(b0 = rep(1, k * nv + n))))
            formula <- y ~ 0 + b0 + f(field, model = spde,
                                      group = field.group, control.group = ctr.g)
        }else{
            stack <- inla.stack(
                data = list(y = y.pp, e = expected),
                A = list(A.pp, 1), 
                effects = list(i = 1:nv, b0 = rep(1, nv + n)))
            formula <- y ~ 0 + b0 + f(i, model = spde)
        }
    }
    ##call to inla
    result <- inla(as.formula(formula), family = "poisson",
                   data = inla.stack.data(stack),
                   E = inla.stack.data(stack)$e,
                   control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
                   control.inla = control.inla,
                   control.fixed = control.fixed,
                   verbose = verbose,
                   ...)
    result
}

#' Function copied from \url{http://www.r-inla.org/spde-book}
#' which is suppied alongside the
#' \href{SPDE gitbook}{https://becarioprecario.bitbucket.io/spde-gitbook/}
#' @source \url{http://www.r-inla.org/spde-book}
inla.mesh.dual <- function(mesh) {
    if (mesh$manifold=='R2') {
        ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
            colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
        library(parallel)
        pls <- mclapply(1:mesh$n, function(i) {
            p <- unique(Reduce('rbind', lapply(1:3, function(k) {
                j <- which(mesh$graph$tv[,k]==i)
                if (length(j)>0) 
                    return(rbind(ce[j, , drop=FALSE],
                                 cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                                       mesh$loc[mesh$graph$tv[j, k], 2] +
                                       mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
                else return(ce[j, , drop=FALSE])
            })))
            j1 <- which(mesh$segm$bnd$idx[,1]==i)
            j2 <- which(mesh$segm$bnd$idx[,2]==i)
            if ((length(j1)>0) | (length(j2)>0)) {
                p <- unique(rbind(mesh$loc[i, 1:2], p,
                                  mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                                  mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                                  mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
                yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
                xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
            }
            else {
                yy <- p[,2]-mesh$loc[i, 2]
                xx <- p[,1]-mesh$loc[i, 1]
            }
            sp::Polygon(p[order(atan2(yy,xx)), ])
        })
        return(sp::SpatialPolygons(lapply(1:mesh$n, function(i)
            sp::Polygons(list(pls[[i]]), i))))
    }
    else stop("It only works for R2!")
}
