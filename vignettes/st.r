## ---- libs
library(stelfi)

## ---- data_iq
data(terrorism)

## ---- maps
iq <- maps::map("world","iraq",fill = TRUE,plot = FALSE)
iq_sp <- maptools::map2SpatialPolygons(iq, IDs = "Iraq",
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))
## transform to utm
iq_sp <-  spTransform(iq_sp,
                      CRS("+proj=utm +zone=38 +datum=WGS84"))
sp <- terrorism
coordinates(sp) <- c("longitude","latitude")
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84")
sp <-  spTransform(sp,CRS("+proj=utm +zone=38 +datum=WGS84"))

## ---- plot
plot(iq_sp)
plot(sp, add = TRUE, col = "grey")
baghdad <- data.frame(y = 33.3152, x = 44.3661)
coordinates(baghdad) <- c("x","y");proj4string(baghdad) <- CRS("+proj=longlat +datum=WGS84")
baghdad <- spTransform(baghdad,CRS("+proj=utm +zone=38 +datum=WGS84"))
text(baghdad$x, baghdad$y, col = 2, "Baghdad")
mosul <- data.frame(y = 36.3489, x = 43.1577)
coordinates(mosul) <- c("x","y");proj4string(mosul) <- CRS("+proj=longlat +datum=WGS84")
mosul <- spTransform(mosul,CRS("+proj=utm +zone=38 +datum=WGS84"))
text(mosul$x, mosul$y, col = 2, "Mosul")

## ---- mesh
## mesh max.edge on the same scale as the coords
mesh <- inla.mesh.2d(loc.domain = broom::tidy(iq_sp)[,1:2],
                     max.edge = c(50000,75000))
## ----  model_inla
## AR(1) over years, 2013--2017
temp <- (sp$iyear - min(sp$iyear)) + 1
fit_inla <- fit_lgcp_inla(locs = coordinates(sp), sp = iq_sp,
                     mesh = mesh, temp = temp,
                     return.attributes = TRUE)
summary(fit_inla)

## ---- fixed_effects
fit_inla$summary.fixed
## areas/weights at each mesh nodes, only returned if return.attribures = TRUE
weights <- attributes(fit_inla)$weights

## ---- random_effects
fields <- stelfi::get_fields(fit_inla, mesh,t = 5, mean = TRUE)

par(mfrow = c(2,3))
for(i in 1:5){
    rf <- fields[[1]][[i]]
    show_field(rf, mesh, dims = c(300,300),
               col = RColorBrewer::brewer.pal(9, "Blues"), sp = iq_sp,
               rast = FALSE, main = 2012 + i)
}

