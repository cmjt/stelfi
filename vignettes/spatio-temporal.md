To install the `R` package `stelfi` run
`devtools::install_github("cmjt/stelfi")`. Use
`devtools::install_github("cmjt/stelfi",build_vignettes = TRUE)` if you
want to build vignettes at the same time.

    library(stelfi)

A point process model for terrorism attacks perpetrated by ISIL in Iraq, 2013--2017
-----------------------------------------------------------------------------------

    data(terrorism)

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

![Point pattern of attacks by ISIL in Iraq,
2013--2017](spatio-temporal_files/figure-markdown_strict/plot-1.png)

    ## mesh max.edge on the same scale as the coords
    mesh <- inla.mesh.2d(loc.domain = broom::tidy(iq_sp)[,1:2],
                         max.edge = c(50000,75000))

### `INLA`

    ## AR(1) over years, 2013--2017
    temp <- (sp$iyear - min(sp$iyear)) + 1
    fit <- fit_lgcp_inla(mesh = mesh, locs = coordinates(sp), sp = iq_sp,
                         temp = temp,
                         return.attributes = TRUE)
    summary(fit)

    ## 
    ## Call:
    ##    c("inla(formula = as.formula(formula), family = \"poisson\", data = 
    ##    inla.stack.data(stack), ", " E = inla.stack.data(stack)$e, verbose = 
    ##    verbose, control.predictor = list(A = inla.stack.A(stack), ", " compute 
    ##    = TRUE), control.inla = control.inla, control.fixed = control.fixed)" ) 
    ## Time used:
    ##     Pre = 10.5, Running = 95.8, Post = 25.1, Total = 131 
    ## Fixed effects:
    ##       mean    sd 0.025quant 0.5quant 0.975quant    mode kld
    ## b0 -19.846 0.016    -19.876  -19.846    -19.815 -19.846   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     field SPDE2 model
    ## 
    ## Model hyperparameters:
    ##                      mean    sd 0.025quant 0.5quant 0.975quant   mode
    ## Range for field    12.254 1.356      9.786   12.189     15.099 12.071
    ## Stdev for field     0.394 0.037      0.327    0.392      0.472  0.387
    ## GroupRho for field  0.779 0.023      0.727    0.781      0.816  0.788
    ## 
    ## Expected number of effective parameters(stdev): 1.03(0.00)
    ## Number of equivalent replicates : 9576.83 
    ## 
    ## Marginal log-Likelihood:  -87738.36 
    ## Posterior marginals for the linear predictor and
    ##  the fitted values are computed

    fit$summary.fixed

    ##         mean         sd 0.025quant  0.5quant 0.975quant      mode          kld
    ## b0 -19.84568 0.01561374  -19.87633 -19.84568  -19.81505 -19.84568 1.970312e-21

    ## areas/weights at each mesh nodes, only returned if return.attribures = TRUE
    weights <- attributes(fit)$weights

    fields <- stelfi::get_fields(fit, mesh,t = 5, mean = TRUE)

    par(mfrow = c(2,3))
    for(i in 1:5){
        rf <- fields[[1]][[i]]
        show_field(rf, mesh, dims = c(300,300),
                   col = RColorBrewer::brewer.pal(9, "Blues"), sp = iq_sp,
                   rast = FALSE, main = 2012 + i)
    }

![Estimated mean of the random fields,
2013--2017](spatio-temporal_files/figure-markdown_strict/random_effects-1.png)

### `TMB`

**TODO**
