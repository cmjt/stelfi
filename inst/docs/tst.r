library(stelfi)
murders_sp <- subset(murders_nz, murders_nz$Region == "Waikato")
## project longitude & latitude to NZTMs
coordinates(murders_sp) <- c("Longitude","Latitude")
proj4string(murders_sp) <- CRS("+proj=longlat +datum=WGS84")
murders_sp <-  spTransform(murders_sp, 
                           CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
waikato <- nz[nz$NAME_1 == "Waikato",]
mesh <- inla.mesh.2d(loc.domain = broom::tidy(waikato)[,1:2],
                     max.edge = c(11000,20000), cutoff = 15000)
## covariate
file <- list.files("data",pattern = ".shp", full = TRUE)
layer <- rgdal::ogrListLayers(file)
pop <- rgdal::readOGR(file, layer = layer)
pop <- spTransform(pop, CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
pop_mesh <- sp::over(SpatialPoints(mesh$loc[,1:2], proj4string = CRS(proj4string(murders_sp))),pop)$Population
pop_obs <- sp::over(murders_sp,pop)$Population
## population density covariate c at mesh nodes and then obs locations 
covs <- data.frame(pop = c(pop_mesh, pop_obs))

fit <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp),covariates = covs, sp = waikato,
                     return.attributes = TRUE)
weights <- attributes(fit)$weights
ins <- which(weights != 0)
en <- exp(as.numeric(fit$summary.fixed[1,1]) +
          as.numeric(fit$summary.fixed[2,1])*covs$pop[1:mesh$n][ins])*weights[ins]
sum(en) ## expected number

fields <- stelfi::get_fields(fit, mesh, mean = TRUE)
grfs <- fields[[1]]
show_field(grfs, mesh, dims = c(300,300),
		 col = RColorBrewer::brewer.pal(9, "Blues"), sp = waikato,
           	 rast = FALSE)

