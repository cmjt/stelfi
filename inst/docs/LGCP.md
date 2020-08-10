To install the `R` package `stelfi` run `devtools::install_github("cmjt/stelfi")`.

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
| Hawke's Bay       |      44|
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

Rate of murders per km<sup>2</sup> across NZ is calculated as 0.0036; there are roughly 60.438 per year. The spatio-temporal murder rate across NZ 2004--2019 is 2.2510^{-4} (rate per km<sup>2</sup> per year).

### Subset to Auckland & transform to NZTM

Transform `data.frame` to `SpatialPointsDataFrame`

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

![](LGCP_files/figure-markdown_github/plot-1.png) *Locations of recorded (n = 108) murders in Waikato 2004--2019*

log-Gaussian Cox process
------------------------

### `stelfi` as a wrapper for `INLA`

#### Creating the mesh

Typically when analysing point pattern data the point locations are not specified as the mesh nodes (i.e., locations are not given as an argument to `inla.mesh.2d()`). Instad we can supply the coordinates of the point pattern window (domain).

``` r
## mesh max.edge on the same scale as the coords (NZTMs)
mesh <- inla.mesh.2d(loc.domain = broom::tidy(waikato)[,1:2],
                     max.edge = c(9000,15000), cutoff = 9000)
```

The SPDE approach for point pattern analysis defines the model at the nodes of the mesh. To fit the log-Cox point process model these points are considered as integration points. The method in Simpson et al. (2016) defines the expected number of events to be proportional to the area around the node (the areas of the polygons in the dual mesh, see below). This means that at the nodes of the mesh with larger triangles, there are also larger expected values.

``` r
weights <- stelfi::get_weights(mesh = mesh, sp = waikato, plot = TRUE)
```

![](LGCP_files/figure-markdown_github/dual%20mesh-1.png)

``` r
## these weights are the areas of each mesh triangle,
## required for the "exposure" aspect of the LGCP.
## the weights are set to zero if outside the study region
## as here the they should have no contribution. 
## the sum of the weights is the area of the study region
```

*Delauney triangulation of the domain (white) overlain on the Voronoi diagram representing the weights (area surrounding) of each mesh node (diamonds). Observations are plotted as circles, mesh nodes outwith the domain are shown in white.*

![](LGCP_files/figure-markdown_github/plot%20weights-1.png) *Voronoi diagram of the weights (shown as areas in km2 around each mesh node).*

``` r
## Spatial only
fit <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp), sp = waikato)
summary(fit)
```

    ## 
    ## Call:
    ##    c("inla(formula = as.formula(formula), family = \"poisson\", data = 
    ##    inla.stack.data(stack), ", " E = inla.stack.data(stack)$e, verbose = 
    ##    verbose, control.predictor = list(A = inla.stack.A(stack), ", " compute 
    ##    = TRUE), control.inla = control.inla, control.fixed = control.fixed)" ) 
    ## Time used:
    ##     Pre = 4.94, Running = 14.2, Post = 0.426, Total = 19.5 
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
    ## Range for field 1.251 1.940      0.116    0.688      5.872 0.281
    ## Stdev for field 0.234 0.269      0.020    0.153      0.939 0.055
    ## 
    ## Expected number of effective parameters(stdev): 1.00(0.00)
    ## Number of equivalent replicates : 1468.99 
    ## 
    ## Marginal log-Likelihood:  -2192.33 
    ## Posterior marginals for the linear predictor and
    ##  the fitted values are computed

``` r
## fixed effects
fit$summary.fixed
```

    ##         mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
    ## b0 -19.23967 0.09624808  -19.42864 -19.23968  -19.05086 -19.23967 6.208903e-25

``` r
## expected number of murders at each mesh node
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
             rast = FALSE, legend = TRUE, FALSE)
```

![](LGCP_files/figure-markdown_github/random%20fields-1.png) *Estimated mean of the assumed Gaussian Markov Random Field*

``` r
## Spatiotemporal
temp <- murders_sp$Year - min(murders_sp$Year) + 1
fit_temp <- fit_lgcp_inla(mesh = mesh, locs = coordinates(murders_sp), sp = waikato,
     temp = temp)
```

### "Raw" `INLA`

**Steps below closely follow this [INLA-SPDE tutorial](https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html)**.

#### Spatial only LGCP

``` r
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
  effects = list(list(b0 = rep(1, nodes + nrow(murders_sp))), 
                 list(i = 1:nodes)),
  tag = 'pp')
## fit model
pp.res <- inla(y ~ 0 + b0 + f(i, model = spde), 
  family = 'poisson', data = inla.stack.data(stk.pp), 
  control.predictor = list(A = inla.stack.A(stk.pp)), 
  E = inla.stack.data(stk.pp)$e)
## fixed effects
pp.res$summary.fixed
```

    ##         mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
    ## b0 -19.23973 0.09621838  -19.43379 -19.23792  -19.05585 -19.23433 2.156986e-07

``` r
## expected number of murders at each mesh node
ins <- which(weights != 0)
en <- exp(as.numeric(pp.res$summary.fixed[1]))*weights[ins]
sum(en) ## expected number across Waikato, observed 108
```

    ## [1] 108.0135

![](LGCP_files/figure-markdown_github/resp-1.png) *Estimated mean of the assumed Gaussian Markov Random Field*

##### Adding a covariate

``` r
## include covariates
## The covariate shapefile used can be downloaded from
## https://koordinates.com/from/datafinder.stats.govt.nz/layer/8437/data/
## the code below assumes a single .shp (above) file is
## in a directory data/ relative to your working directory
file <- list.files("data",pattern = ".shp", full = TRUE)
layer <- rgdal::ogrListLayers(file)
pop <- rgdal::readOGR(file, layer = layer)
```

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "/home/charlotte/Git/stelfi/inst/docs/data/population-by-meshblock-2013-census.shp", layer: "population-by-meshblock-2013-census"
    ## with 46621 features
    ## It has 4 fields

``` r
pop <- spTransform(pop, CRS("+proj=nzmg +lat_0=-41.0 +lon_0=173.0 +x_0=2510000.0 +y_0=6023150.0 +ellps=intl +units=m"))
pop_mesh <- sp::over(SpatialPoints(mesh$loc[,1:2], proj4string = CRS(proj4string(murders_sp))),pop)$Population
## will obviously be NA at mesh nodes outside NZ
```

![](LGCP_files/figure-markdown_github/inference-1.png) *Voronoi diagram of the covariate (population) per mesh node.*

``` r
pop_obs <- sp::over(murders_sp,pop)$Population
## population density covariate c at mesh nodes and then obs locations 
covs <- data.frame(pop = c(pop_mesh, pop_obs))
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

    ##              mean           sd    0.025quant      0.5quant    0.975quant
    ## b0  -19.686379656 0.1265888536 -19.940116344 -19.684573519 -19.442881792
    ## pop   0.005619561 0.0007521601   0.004053218   0.005650749   0.007009032
    ##              mode          kld
    ## b0  -19.680987262 9.994579e-09
    ## pop   0.005713182 2.379353e-06

``` r
## expected number of murders at each mesh node
ins <- which(weights != 0)
en <- exp(as.numeric(pp.cov$summary.fixed[1,1]) +
          as.numeric(pp.cov$summary.fixed[2,1])*covs$pop[1:mesh$n][ins])*weights[ins]
sum(en) ## expected number 
```

    ## [1] 108.5103

#### Spatio-tempoal LGCP

``` r
## space time SPDE
## A set of knots over time needs to be defined in order to fit a
## SPDE spatio-temporal model. It is then used to build a temporal mesh, as follows:

## spatio temporal
k <- length(table(murders_sp$Year))
temp <- murders_sp$Year - min(murders_sp$Year) + 1
## tknots <- seq(min(data$Year), max(data$Year), length = k) ## don't need this as year already discrete
mesh.t <- inla.mesh.1d(1:k)
## spatial spde
spde <- inla.spde2.pcmatern(mesh = mesh,
  prior.range = c(5, 0.01), # P(practic.range < 5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
m <- spde$n.spde
## spatio-temporal projection matrix
Ast <- inla.spde.make.A(mesh = mesh, loc = coordinates(murders_sp),
                        n.group = length(mesh.t$n), group = temp,
                        group.mesh = mesh.t)
## index set
idx <- inla.spde.make.index('s', spde$n.spde, n.group = mesh.t$n)
## spatio-temporal volume
st.vol <- rep(weights, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), m)
y <- rep(0:1, c(k * m, nrow(murders_sp)))
expected <- c(st.vol, rep(0, nrow(murders_sp)))
stk <- inla.stack(
    data = list(y = y, expect = expected), 
    A = list(rbind(Diagonal(n = k * m), Ast), 1), 
    effects = list(idx, list(a0 = rep(1, k * m + nrow(murders_sp)))))
## formula
pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))
form <- y ~ 0 + a0 + f(s, model = spde, group = s.group, 
                       control.group = list(model = 'ar1',
                                            hyper = list(theta = pcrho)))
res <- inla(form, family = 'poisson', 
            data = inla.stack.data(stk), E = expect,
            control.predictor = list(A = inla.stack.A(stk)),
            control.inla = list(strategy = 'adaptive'))
```

### Using `TMB` TODO

``` r
## ensure you've run compile.stelfi()
```
