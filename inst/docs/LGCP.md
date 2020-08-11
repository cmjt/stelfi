To install the `R` package `stelfi` run
`devtools::install_github("cmjt/stelfi")`.

``` r
library(stelfi)
```

NZ murders
==========

The Data
--------

``` r
data(murders_nz)
dim(murders_nz)
```

    ## [1] 967  13

``` r
head(murders_nz)
```

    ##    Latitude Longitude    Sex Age   Date Year                      Cause  Killer
    ## 1 -43.63394  171.6442   Male  41  Jan 5 2004                   stabbing  friend
    ## 2 -43.28563  172.1305   Male  46  Jan 8 2004            pick axe wounds  friend
    ## 3 -36.92575  174.8498   Male   0 Jan 15 2004 asphyxiation (suffocation)  mother
    ## 4 -43.55006  172.6327 Female  46  Feb 1 2004         blunt force trauma partner
    ## 5 -40.73297  175.1195   Male  10  Feb 2 2004                   stabbing  father
    ## 6 -40.73273  175.1193 Female   2  Feb 2 2004                   stabbing  father
    ##                      Name  Full_date    Month          Cause_cat     Region
    ## 1          Donald Linwood 2004-01-05  January     Violent weapon Canterbury
    ## 2             James Weeks 2004-01-08  January     Violent weapon Canterbury
    ## 3 Gabriel Harrison-Taylor 2004-01-15  January           Asphyxia   Auckland
    ## 4   Odette Lloyd-Rangiuia 2004-02-01 February Blunt force trauma Canterbury
    ## 5         Te Hau OCarroll 2004-02-02 February     Violent weapon Wellington
    ## 6        Ngamata OCarroll 2004-02-02 February     Violent weapon Wellington

| Region            |  Number|
|:------------------|-------:|
| Auckland          |     267|
| Bay of Plenty     |      94|
| Canterbury        |     116|
| Gisborne          |      17|
| Hawke’s Bay       |      44|
| Manawatu-Wanganui |      65|
| Marlborough       |       6|
| Nelson            |      16|
| Northland         |      49|
| Otago             |      30|
| Southland         |      19|
| Taranaki          |      27|
| Tasman            |       6|
| Waikato           |     108|
| Wellington        |      89|
| West Coast        |      11|

*Number of murders by province.*

``` r
data(nz) ## SpatialPolygonsDataFrame of NZ (NZTM projection)
area_nz <- sum(raster::area(nz)) ## (m)
area_nzkm2 <- area_nz/1000^2 ## according to Google NZ is 268,021 km2
area_nzkm2
```

    ## [1] 268856.4

``` r
spatial_murder_rate <- nrow(murders_nz)/area_nzkm2 
temporal_murder_rate <- nrow(murders_nz)/length(table(murders_nz$Year)) 
st_murder_rate <- (nrow(murders_nz)/area_nzkm2)/length(table(murders_nz$Year)) 
```

There are on average 60.4 murders per year. The rate of murders per
km<sup>2</sup> across NZ is roughly 0.004.

### Subset data to Waikato & transform to NZTM

Transform `data.frame` to `SpatialPointsDataFrame` and project Longitude
and Latitude to NZTM.

``` r
murders_sp <- subset(murders_nz, murders_nz$Region == "Waikato")
## project longitude & latitude to NZTMs
coordinates(murders_sp) <- c("Longitude","Latitude")
proj4string(murders_sp) <- CRS("+proj=longlat +datum=WGS84")
murders_sp <-  spTransform(murders_sp, 
                           CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
waikato <- nz[nz$NAME_1 == "Waikato",]
```

| Region             |  Number|
|:-------------------|-------:|
| Asphyxia           |       2|
| Blunt force trauma |      35|
| Car crash          |      13|
| Drowning           |       1|
| Fire               |       4|
| Other              |      13|
| Violent weapon     |      40|

*Number of murders by category in Waikato 2004–2019.*

![](LGCP_files/figure-markdown_github/plot-1.png) *Locations of recorded
(n = 108) murders in Waikato 2004–2019*

log-Gaussian Cox process
------------------------

### `stelfi` as a wrapper for `INLA`

#### Creating the mesh

Typically when analysing point pattern data the point locations are not
specified as the mesh nodes (i.e., locations are not given as an
argument to `inla.mesh.2d()`). Instad we can supply the coordinates of
the point pattern window (domain).

``` r
## mesh max.edge on the same scale as the coords (NZTMs)
mesh <- inla.mesh.2d(loc.domain = broom::tidy(waikato)[,1:2],
                     max.edge = c(11000,20000), cutoff = 15000)
```

#### Model fitting

``` r
## Spatial only
fit <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp), sp = waikato,
                     return.attributes = TRUE)
summary(fit)
```

    ## 
    ## Call:
    ##    c("inla(formula = as.formula(formula), family = \"poisson\", data = 
    ##    inla.stack.data(stack), ", " E = inla.stack.data(stack)$e, verbose = 
    ##    verbose, control.predictor = list(A = inla.stack.A(stack), ", " compute 
    ##    = TRUE), control.inla = control.inla, control.fixed = control.fixed)" ) 
    ## Time used:
    ##     Pre = 0.542, Running = 1.28, Post = 0.167, Total = 1.99 
    ## Fixed effects:
    ##      mean    sd 0.025quant 0.5quant 0.975quant   mode kld
    ## b0 -19.24 0.096    -19.429   -19.24    -19.051 -19.24   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     field SPDE2 model
    ## 
    ## Model hyperparameters:
    ##                  mean    sd 0.025quant 0.5quant 0.975quant  mode
    ## Range for field 0.886 1.063      0.100    0.567       3.62 0.255
    ## Stdev for field 0.272 0.294      0.029    0.185       1.05 0.078
    ## 
    ## Expected number of effective parameters(stdev): 1.00(0.00)
    ## Number of equivalent replicates : 876.40 
    ## 
    ## Marginal log-Likelihood:  -2192.33 
    ## Posterior marginals for the linear predictor and
    ##  the fitted values are computed

``` r
## fixed effects
fit$summary.fixed
```

    ##         mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
    ## b0 -19.23967 0.09624817  -19.42864 -19.23968  -19.05086 -19.23967 6.982366e-26

``` r
## expected number of murders at each mesh node
## areas/weights at each mesh nodes, only returned if return.attribures = TRUE
weights <- attributes(fit)$weights
ins <- which(weights != 0)
en <- exp(as.numeric(fit$summary.fixed[1]))*weights[ins]
sum(en) ## expected number across Waikato, observed 108
```

    ## [1] 108.0193

``` r
fields <- stelfi::get_fields(fit, mesh, mean = TRUE)
grfs <- fields[[1]]
show_field(grfs, mesh, dims = c(300,300),
         col = RColorBrewer::brewer.pal(9, "Blues"), sp = waikato,
             rast = FALSE, legend = TRUE, legend.only = FALSE)
```

![](LGCP_files/figure-markdown_github/random%20fields-1.png) *Estimated
mean of the assumed Gaussian Markov Random Field*

``` r
fields <- stelfi::get_fields(fit, mesh, mean = FALSE)
grfsd <- fields[[1]]
show_field(grfsd, mesh, dims = c(300,300),
         col = RColorBrewer::brewer.pal(9, "Blues"), sp = waikato,
             rast = FALSE)
```

![](LGCP_files/figure-markdown_github/random%20fields%20sd-1.png)
*Standard error of the assumed Gaussian Markov Random Field*

##### Including a covariate

``` r
## include covariates
## The covariate shapefile used can be downloaded from
## https://koordinates.com/from/datafinder.stats.govt.nz/layer/8437/data/
## the code below assumes a single .shp (above) file is
## in a directory data/ relative to your working directory
file <- list.files("data",pattern = ".shp", full = TRUE)
layer <- rgdal::ogrListLayers(file)
pop <- rgdal::readOGR(file, layer = layer)
pop <- spTransform(pop, CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
pop_mesh <- sp::over(SpatialPoints(mesh$loc[,1:2], proj4string = CRS(proj4string(murders_sp))),pop)$Population
## will obviously be NA at mesh nodes outside NZ
pop_obs <- sp::over(murders_sp,pop)$Population
## population density covariate c at mesh nodes and then obs locations 
covs <- data.frame(pop = c(pop_mesh, pop_obs))
```

``` r
fit <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp),covariates = covs, sp = waikato,
                     return.attributes = TRUE)
weights <- attributes(fit)$weights
ins <- which(weights != 0)
en <- exp(as.numeric(fit$summary.fixed[1,1]) +
          as.numeric(fit$summary.fixed[2,1])*covs$pop[1:mesh$n][ins])*weights[ins]
sum(en) ## expected number 
```

    ## [1] 108.021

``` r
fields <- stelfi::get_fields(fit, mesh, mean = TRUE)
grfs <- fields[[1]]
show_field(grfs, mesh, dims = c(300,300),
         col = RColorBrewer::brewer.pal(9, "Blues"), sp = waikato,
             rast = FALSE)
```

![](LGCP_files/figure-markdown_github/random%20fields%20cov-1.png)
*Estimated mean of the assumed Gaussian Markov Random Field*

### “Raw” `INLA`

**Steps below closely follow this [INLA-SPDE
tutorial](https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html)**.

The SPDE approach for point pattern analysis defines the model at the
nodes of the mesh. To fit the log-Cox point process model these points
are considered as integration points. The method in Simpson et al.
(2016) defines the expected number of events to be proportional to the
area around the node (the areas of the polygons in the dual mesh, see
below). This means that at the nodes of the mesh with larger triangles,
there are also larger expected values.

``` r
## these weights are the areas of each mesh triangle,
## required for the "exposure" aspect of the LGCP.
## the weights are set to zero if outside the study region
## as here the they should have no contribution. 
## the sum of the weights is the area of the study region
weights <- stelfi::get_weights(mesh = mesh, sp = waikato, plot = TRUE)
```

![](LGCP_files/figure-markdown_github/dual%20mesh-1.png) *Delauney
triangulation of the domain (white) overlain on the Voronoi diagram
representing the weights (area surrounding) of each mesh node
(diamonds). Observations are plotted as circles, mesh nodes outwith the
domain are shown in white.*

![](LGCP_files/figure-markdown_github/plot%20weights-1.png) *Voronoi
diagram of the weights (shown as areas in km2 around each mesh node).*

``` r
## number of mesh nodes
nodes <- mesh$n
## define model
spde <- inla.spde2.pcmatern(mesh = mesh,
  # PC-prior on range: P(practic.range < 5) = 0.9
  prior.range = c(5, 0.9),
  # PC-prior on sigma: P(sigma > 1) = 0.005
  prior.sigma = c(1, 0.005))
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
  effects = list(list(b0 = rep(1, nodes + nrow(murders_sp))), 
                 list(i = 1:nodes)),
  tag = 'pp')
## fit model
pp.res <- inla(y ~ 0 + b0 + f(i, model = spde), 
  family = 'poisson', data = inla.stack.data(stk.pp), 
  control.predictor = list(A = inla.stack.A(stk.pp)), 
  E = inla.stack.data(stk.pp)$e)
summary(pp.res)
```

    ## 
    ## Call:
    ##    c("inla(formula = y ~ 0 + b0 + f(i, model = spde), family = 
    ##    \"poisson\", ", " data = inla.stack.data(stk.pp), E = 
    ##    inla.stack.data(stk.pp)$e, ", " control.predictor = list(A = 
    ##    inla.stack.A(stk.pp)))") 
    ## Time used:
    ##     Pre = 1.51, Running = 2.72, Post = 0.153, Total = 4.39 
    ## Fixed effects:
    ##      mean    sd 0.025quant 0.5quant 0.975quant    mode kld
    ## b0 -19.24 0.096    -19.434  -19.238    -19.056 -19.234   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     i SPDE2 model
    ## 
    ## Model hyperparameters:
    ##              mean    sd 0.025quant 0.5quant 0.975quant  mode
    ## Range for i 0.886 1.063      0.100    0.567       3.62 0.255
    ## Stdev for i 0.272 0.294      0.029    0.185       1.05 0.078
    ## 
    ## Expected number of effective parameters(stdev): 1.00(0.00)
    ## Number of equivalent replicates : 876.40 
    ## 
    ## Marginal log-Likelihood:  -2192.33

``` r
## fixed effects
pp.res$summary.fixed
```

    ##         mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
    ## b0 -19.23973 0.09621868  -19.43378 -19.23792  -19.05584 -19.23435 2.125873e-07

``` r
## expected number of murders at each mesh node
ins <- which(weights != 0)
en <- exp(as.numeric(pp.res$summary.fixed[1]))*weights[ins]
sum(en) ## expected number across Waikato, observed 108
```

    ## [1] 108.0136

![](LGCP_files/figure-markdown_github/resp-1.png) *Estimated mean of the
assumed Gaussian Markov Random Field*

##### Adding a covariate

``` r
## data stack
stk.cov <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nodes + nrow(murders_sp)), pop = covs$pop), 
                 list(i = 1:nodes)),
  tag = 'pp')
## fit model
pp.cov <- inla(y ~ 0 + b0 + pop + f(i, model = spde), 
  family = 'poisson', data = inla.stack.data(stk.cov), 
  control.predictor = list(A = inla.stack.A(stk.cov)), 
  E = inla.stack.data(stk.pp)$e)
## coefficients of the fixed effects
pp.cov$summary.fixed
```

    ##             mean          sd    0.025quant     0.5quant   0.975quant
    ## b0  -20.20124163 0.181506672 -20.569372105 -20.19706941 -19.85652196
    ## pop   0.01139917 0.001410218   0.008623026   0.01140167   0.01415873
    ##             mode          kld
    ## b0  -20.18879075 5.628601e-07
    ## pop   0.01140677 1.027338e-06

``` r
## expected number of murders at each mesh node
ins <- which(weights != 0)
en <- exp(as.numeric(pp.cov$summary.fixed[1,1]) +
          as.numeric(pp.cov$summary.fixed[2,1])*covs$pop[1:mesh$n][ins])*weights[ins]
sum(en) ## expected number 
```

    ## [1] 108.2688

### Using `TMB` TODO

``` r
## ensure you've run compile.stelfi()
```
