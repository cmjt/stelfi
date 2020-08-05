library(stelfi)
murders_sp <- murders_nz
## project longitude & latitude to NZTMs
coordinates(murders_sp) <- c("Longitude","Latitude")
proj4string(murders_sp) <- CRS("+proj=longlat +datum=WGS84")
murders_sp <-  spTransform(murders_sp, 
                           CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
mesh <- inla.mesh.2d(loc.domain = coordinates(nz) ,
                     max.edge = c(86000, 100000), cutoff = 5000)
##########################################
## covariate
## include covariates
## let's assume you have a covariate shapefile
## shape file from https://koordinates.com/layer/7322-new-zealand-population-density-by-meshblock/
file <- list.files("data",pattern = ".shp", full = TRUE)
layer <- rgdal::ogrListLayers(file)
pop <- rgdal::readOGR(file, layer = layer)
pop <- spTransform(pop, CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
pop_mesh <- sp::over(SpatialPoints(mesh$loc[,1:2], proj4string = CRS(proj4string(murders_sp))),pop)
## will obviously be NA at mesh nodes outside NZ, don't worry
pop_obs <- sp::over(murders_sp,pop)
## population density covariate c at mesh nodes and then obs locations 
covs <- data.frame(pop = c(pop_mesh$pop_densit, pop_obs$pop_densit))
##########################################
weights <- stelfi:::get_weights(mesh, nz, TRUE)
## number of mesh nodes
nodes <- mesh$n
## define model
spde <- inla.spde2.pcmatern(mesh = mesh,
  # PC-prior on range: P(practic.range < 0.05) = 0.01
  prior.range = c(0.05, 0.01),
  # PC-prior on sigma: P(sigma > 1) = 0.01
  prior.sigma = c(1, 0.01))
## vector for observations
y.pp <- rep(0:1, c(nodes, nrow(murders_sp)))
## exposure (E)
e.pp <- c(weights, rep(0, nrow(murders_sp)))
## integration points
imat <- Diagonal(nodes, rep(1, nodes))
## projection matrix for observed points
lmat <- inla.spde.make.A(mesh, coordinates(murders_sp))
## entire projection matrix
A.pp <- rbind(imat, lmat)
## data stack
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nodes + nrow(murders_sp)), pop = covs$pop), 
                 list(i = 1:nodes)),
  tag = 'pp')
## fit model
pp.res <- inla(y ~ 0 + b0 + pop + f(i, model = spde), 
  family = 'poisson', data = inla.stack.data(stk.pp), 
  control.predictor = list(A = inla.stack.A(stk.pp)), 
  E = inla.stack.data(stk.pp)$e)
## ## fixed effects
## pp.res$summary.fixed
## ## expected number of murders at each mesh node
## ins <- which(weights != 0)
## en <- exp(as.numeric(pp.res$summary.fixed[1]))*weights[ins]
## sum(en) ## expected number across NZ, observed 967
## dual_mesh <- stelfi:::inla.mesh.dual(mesh)


## ## spatio temporal
## k <- length(table(murders_sp$Year))
## temp <- murders_sp$Year - min(murders_sp$Year) + 1
## ## tknots <- seq(min(data$Year), max(data$Year), length = k) ## don't need this as year already discrete
## mesh.t <- inla.mesh.1d(1:k)
## ## spatial spde
## spde <- inla.spde2.pcmatern(mesh = mesh,
##   prior.range = c(5, 0.01), # P(practic.range < 5) = 0.01
##   prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
## m <- spde$n.spde
## ## spatio-temporal projection matrix
## Ast <- inla.spde.make.A(mesh = mesh, loc = coordinates(murders_sp),
##                         n.group = length(mesh.t$n), group = temp,
##                         group.mesh = mesh.t)
## ## index set
## idx <- inla.spde.make.index('s', spde$n.spde, n.group = mesh.t$n)
## ## spatio-temporal volume
## st.vol <- rep(weights, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), m)
## y <- rep(0:1, c(k * m, nrow(murders_sp)))
## expected <- c(st.vol, rep(0, nrow(murders_sp)))
## stk <- inla.stack(
##     data = list(y = y, expect = expected), 
##     A = list(rbind(Diagonal(n = k * m), Ast), 1), 
##     effects = list(idx, list(a0 = rep(1, k * m + nrow(murders_sp)))))
## ## formula
## pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))
## form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
##                        control.group = list(model = 'ar1',
##                                             hyper = list(theta = pcrho)))
## res <- inla(form, family = 'poisson', 
##             data = inla.stack.data(stk), E = expect,
##             control.predictor = list(A = inla.stack.A(stk)),
##             control.inla = list(strategy = 'adaptive'))

## stelfi
st <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp), sp = nz)

## fit model with covariate
scov <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp),covariate =  covs,sp = nz)


temp <- murders_sp$Year - min(murders_sp$Year) + 1
## stemp <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp), sp = nz,
##                        temp = temp)

int <- split(x <- replicate(max(temp), rep(1,mesh$n)),col(x)); names(int) <- rep("intercept",length(int))
## fit <- fit_lgcp_tmb(locs = coordinates(murders_sp), temp.idx = as.factor(temp),
##                     mesh = mesh, parameters = list(beta = 1, log_kappa = -0.05, rho = 0.6),
##                     covs = int,
##                     sp = nz)

##****************### need to make mesh with locs
mesh <- inla.mesh.2d(loc = coordinates(murders_sp),loc.domain = coordinates(nz) ,
                     max.edge = c(86000, 100000), cutoff = 5000)
resp <- matrix(NA,ncol = length(table(temp)), nrow = mesh$n)
temp.idx <- rep(1:2, times = c(500,467)) #as.factor(temp)
w.loc <- split(mesh$idx$loc,temp.idx)
w.c <- lapply(w.loc,table)
for(i in 1:length(w.c)){
    resp[,i] <- numeric(mesh$n)
    count <- as.vector(w.c[[i]])
    resp[,i][unique(w.loc[[i]])] <- count
}
covs <- lapply(int, as.matrix)
data <- list(y = resp, covariates = covs, tsteps = length(table(temp)))
spde <- inla.spde2.matern(mesh = mesh,alpha = 2)
data$spde <- spde$param.inla[c("M0","M1","M2")]
w <- get_weights(mesh = mesh,sp = nz,FALSE)
data$area <- w#*c(Matrix::diag(data$spde$M0))

## params
parameters <- list(beta = c(2), log_rho = -0.1, log_sigma = -0.02 , log_kappa = -1.5,
                   field =  as.matrix(matrix(0,nrow = mesh$n, ncol = length(table(temp.idx)))))

## compile.stelfi()
## dll.stelfi()
fit <- TMB::MakeADFun(data, parameters, DLL = "lgcpar1", random = c("field"))
## Fitting the model.
fit.fixed <- optim(fit$par, fit$fn, fit$gr)
## Getting sdreport.
sdrep.fixed <- sdreport(obj.fixed)
