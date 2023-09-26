test_that("meshmetrics()", {
    data(horse_mesh, package = "stelfi")
    metrics <- meshmetrics(horse_mesh)
    expect_equal(metrics$quality[4],
                 0.8843856,
                 tolerance = 0.01)
})
test_that("Simple Hawkes model diagnostics", {
    ## NIWA retweets
    data(retweets_niwa, package = "stelfi")
    times <- unique(sort(as.numeric(difftime(retweets_niwa,
                                             min(retweets_niwa),
                                             units = "mins"))))
    params <- c(mu = 9, alpha = 3, beta = 10)
    fit <- fit_hawkes(times = times, parameters = params)
    c <- show_hawkes_GOF(obj = fit, plot = FALSE, return_values = TRUE)
    expect_equal(c[[1]][1], 0.15848, tolerance = 0.2)
})
test_that("Get LGCP fields (spatial)", {
    skip_on_cran()
    if(requireNamespace("fmesher")){
        data(xyt, package = "stelfi")
        domain <- sf::st_as_sf(xyt$window)
        locs <- data.frame(x = xyt$x, y = xyt$y)
        bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
        smesh <-  fmesher::fm_mesh_2d(boundary = bnd,
                                    max.edge = 0.75, cutoff = 0.3)
        fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
                        parameters = c(beta = 0, log_tau = log(1),
                                       log_kappa = log(1)))
        f <- get_fields(fit, smesh)
        expect_equal(f[1], 0.149, tolerance = 0.2)
    }
})
test_that("Get LGCP fields (spatiotemporal)", {
    skip_on_cran()
    if(requireNamespace("fmesher")){
        data(xyt, package = "stelfi")
        domain <- sf::st_as_sf(xyt$window)
        ndays <- 2
        locs <- data.frame(x = xyt$x, y = xyt$y, t = xyt$t)
        bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
        w0 <- 2
        smesh <- fmesher::fm_mesh_2d(boundary = bnd,
                                    max.edge = 0.75, cutoff = 0.3)
        tmesh <- fmesher::fm_mesh_1d(seq(0, ndays, by = w0))
        fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh, tmesh = tmesh,
                        parameters = c(beta = 0, log_tau = log(1),
                                       log_kappa = log(1), atanh_rho = 0.2))
        f <- get_fields(fit, smesh)
        expect_equal(f[1], 0.4561646, tolerance = 0.5)
    }
})
test_that("stelfi:::segments()", {
    data(horse_mesh, package = "stelfi")
    seg <- segments(horse_mesh)
    expect_equal(seg[1, 5],
                 1.684799,
                 tolerance = 0.01)
})
