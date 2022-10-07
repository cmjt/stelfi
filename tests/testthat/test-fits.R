test_that("Hawkes simulation", {
    ## sim_hawkes example
    set.seed(1234)
    times <- sim_hawkes(mu = 0.3, alpha = 4, beta = 5)
    expect_equal(round(times[1:6], 2),
                 c(4.15,  4.95,  7.80, 14.64, 18.69, 33.61),
                 tolerance = 0.01)
})
test_that("Simple Hawkes model fitting", {
    ## NIWA retweets
    data(retweets_niwa)
    times <- unique(sort(as.numeric(difftime(retweets_niwa,
                                             min(retweets_niwa),
                                             units = "mins"))))
    params <- c(mu = 9, alpha = 3, beta = 10)
    fit <- fit_hawkes(times = times, parameters = params)
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 0.06328, tolerance = 1)
    expect_equal(pars[2], 0.07597, tolerance = 1)
    expect_equal(pars[3], 0.07911, tolerance = 1)
})
test_that("Non-homogeneous Hawkes model fitting", {
    set.seed(1)
    times <- sim_hawkes(mu = 0.3, alpha = 4, beta = 5)
    background <- function(params, times) {
        A <- exp(params[[1]])
        B <- stats::plogis(params[[2]]) * A
        return(A + B * sin(times))
    }
    background_integral <- function(params, x) {
        A <- exp(params[[1]])
        B <- stats::plogis(params[[2]]) * A
        return((A * x) - B * cos(x))
    }
    param <- list(alpha = 0.5, beta = 1.5)
    background_param <- list(1, 1)
    fit <- fit_hawkes_cbf(times = times, parameters = param,
                          background = background,
                          background_integral = background_integral,
                          background_parameters = background_param)
    estA <- exp(fit$background_parameters[1])
    estB <- plogis(fit$background_parameters[2]) * exp(fit$background_parameters[1])
    expect_equal(estA, 0.2387370, tolerance = 0.1)
    expect_equal(estB, 0.02600446, tolerance = 0.1)
})
test_that("LGCP model fitting (spatial)", {
    skip_on_cran()
    stelfi_load_inla()
    data(xyt, package = "stelfi")
    domain <- sf::st_as_sf(xyt$window)
    locs <- data.frame(x = xyt$x, y = xyt$y)
    bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    smesh <- INLA::inla.mesh.2d(boundary = bnd,
                                max.edge = 0.75, cutoff = 0.3)
    fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
                                parameters = c(beta = 0, log_tau = log(1),
                                               log_kappa = log(1)))
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 2.45, tolerance = 0.5)
    expect_equal(pars[2], -1.32, tolerance = 0.1)
    expect_equal(pars[3], 0.95, tolerance = 0.1)
})
test_that("LGCP model fitting (spatiotemporal)", {
    skip_on_cran()
    stelfi_load_inla()
    require(maptools)
    data(xyt, package = "stelfi")
    domain <- sf::st_as_sf(xyt$window)
    ndays <- 2
    locs <- data.frame(x = xyt$x, y = xyt$y, t = xyt$t)
    bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    w0 <- 2
    smesh <- INLA::inla.mesh.2d(boundary = bnd,
                                max.edge = 0.75, cutoff = 0.3)
    tmesh <- INLA::inla.mesh.1d(seq(0, ndays, by = w0))
    fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh, tmesh = tmesh,
                    parameters = c(beta = 0, log_tau = log(1),
                                   log_kappa = log(1), atanh_rho = 0.2))
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 0.31, tolerance = 1)
    expect_equal(pars[2], 1.68, tolerance = 1)
    expect_equal(pars[3], -1.06, tolerance = 1)
})
test_that("LGCP model fitting (marked)", {
    skip_on_cran()
    stelfi_load_inla()
    data(marked, package = "stelfi")
    loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
    domain <-  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(loc.d))))
    smesh <- INLA::inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1),
                                max.edge = c(0.3, 0.7), cutoff = 0.05)
    locs <- cbind(x = marked$x, y = marked$y)
    marks <- cbind(m1 = marked$m1) ## Gaussian
    parameters <- list(betamarks = matrix(0, nrow = 1, ncol = ncol(marks)),
                       log_tau = log(1), log_kappa = log(1),
                       marks_coefs_pp = rep(0, ncol(marks)), betapp = 0)

    fit <- fit_mlgcp(locs = locs, marks = marks,
                     sf = domain, smesh = smesh,
                     parameters = parameters, methods = 0,
                     fields = 0)
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 9.90, tolerance = 0.1)
    expect_equal(pars[2], 2.74, tolerance = 0.1)
    expect_equal(pars[3], -0.279, tolerance = 0.01)
})
test_that("Spatial self-exciting (no GMRF)", {
    skip_on_cran()
    stelfi_load_inla()
    data(xyt, package = "stelfi")
    N <- 50
    locs <- data.frame(x = xyt$x[1:N], y = xyt$y[1:N])
    times <- xyt$t[1:N]
    domain <- sf::st_as_sf(xyt$window)
    bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    smesh <- INLA::inla.mesh.2d(boundary = bnd,
                                max.edge = 0.75, cutoff = 0.3)
    param <- list(mu = 3, alpha = 1, beta = 3, xsigma = 0.2,
                  ysigma = 0.2, rho = 0.8)
    fit <- fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh, parameters = param) 
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 0.3, tolerance = 0.2)
    expect_equal(pars[2], -1.1, tolerance = 0.2)
    expect_equal(pars[7], -0.17, tolerance = 0.2)
})
test_that("Spatial self-exciting (GMRF)", {
    skip_on_cran()
    stelfi_load_inla()
    data(xyt, package = "stelfi")
    N <- 50
    locs <- data.frame(x = xyt$x[1:N], y = xyt$y[1:N])
    times <- xyt$t[1:N]
    domain <- sf::st_as_sf(xyt$window)
    bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    smesh <- INLA::inla.mesh.2d(boundary = bnd,
                                max.edge = 0.75, cutoff = 0.3)
    param <- list( mu = 5, alpha = 1, beta = 3, kappa = 0.9, tau = 1, xsigma = 0.2, ysigma = 0.2, rho = 0.8)
    fit <- fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh, parameters = param, GMRF = TRUE)
    pars <- as.numeric(get_coefs(fit)[, 1])
    expect_equal(pars[1], 0.003, tolerance = 0.02)
    expect_equal(pars[2], -5.88, tolerance = 0.2)
})
