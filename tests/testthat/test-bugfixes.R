## Regression tests for bug fixes in stelfi 1.0.4

## negative log-likelihood of a multivariate Hawkes process (decay rate of the excited stream)
mhawkes_nll_r <- function(times, stream, mu, alpha, beta) {
    s <- as.integer(factor(stream))
    n <- length(times)
    last <- times[n]
    ll <- 0
    for (i in seq_len(n)) {
        j <- s[i]
        prev <- seq_len(i - 1)
        ll <- ll + log(mu[j] + sum(alpha[j, s[prev]] * exp(-beta[j] * (times[i] - times[prev]))))
    }
    for (j in seq_along(mu)) {
        prev <- seq_len(n - 1) ## events before the last event
        ll <- ll - mu[j] * last -
            sum(alpha[j, s[prev]] / beta[j] * (1 - exp(-beta[j] * (last - times[prev]))))
    }
    -ll
}

test_that("Multivariate Hawkes likelihood matches direct calculation", {
    data(multi_hawkes, package = "stelfi")
    mu <- c(0.2, 0.3)
    alpha <- matrix(c(0.5, 0.05, 0.3, 0.4), nrow = 2)
    beta <- c(0.7, 1.2)
    fit <- fit_mhawkes(times = multi_hawkes$times, stream = multi_hawkes$stream,
                       parameters = list(mu = mu, alpha = alpha, beta = beta))
    par <- c(log(mu), qlogis(as.vector(sweep(alpha, 1, beta, "/"))), log(beta))
    expect_equal(fit$fn(par),
                 mhawkes_nll_r(multi_hawkes$times, multi_hawkes$stream, mu, alpha, beta),
                 tolerance = 1e-6)
    ## last event in either stream
    tms <- multi_hawkes$times[1:50]
    str <- multi_hawkes$stream[1:50]
    str[50] <- setdiff(unique(multi_hawkes$stream), str[49])[1]
    fit2 <- fit_mhawkes(times = tms, stream = str,
                        parameters = list(mu = mu, alpha = alpha, beta = beta))
    expect_equal(fit2$fn(par), mhawkes_nll_r(tms, str, mu, alpha, beta), tolerance = 1e-6)
})

test_that("Multivariate Hawkes intensity uses alpha[j, k] and beta[j]", {
    times <- c(1, 2)
    stream <- c("A", "B")
    mu <- c(0.1, 0.2)
    alpha <- matrix(c(0.5, 0.1, 0.3, 0.4), nrow = 2) ## alpha[1, 2] = 0.3, alpha[2, 1] = 0.1
    beta <- c(1, 2)
    lam <- stelfi:::multi_hawkes_intensity(times = times, mu = mu, alpha = alpha, beta = beta,
                                           p = 3, stream = stream)
    expect_equal(lam[[1]], 0.1 + 0.5 * exp(-1 * 2) + 0.3 * exp(-1 * 1))
    expect_equal(lam[[2]], 0.2 + 0.1 * exp(-2 * 2) + 0.4 * exp(-2 * 1))
})

test_that("Spatiotemporal LGCP: ordering, empty intervals, and knot range", {
    skip_on_cran()
    disc <- sf::st_sf(geometry = sf::st_geometry(sf::st_as_sf(spatstat.geom::disc())))
    smesh <- fmesher::fm_mesh_2d(boundary = fmesher::fm_as_segm(disc), max.edge = 0.3, cutoff = 0.1)
    tmesh <- fmesher::fm_mesh_1d(0:3)
    set.seed(1)
    n <- 200
    r <- sqrt(runif(n)) * 0.9; th <- runif(n, 0, 2 * pi)
    locs <- data.frame(x = r * cos(th), y = r * sin(th), t = runif(n, 0, 3))
    locs$t[1] <- 3 ## an event on the last knot
    sv <- c(beta = 3, log_tau = 0, log_kappa = 0, atanh_rho = 0.2)
    srt <- locs[order(locs$t), ]
    fit_sorted <- fit_lgcp(locs = srt, sf = disc, smesh = smesh, tmesh = tmesh, parameters = sv)
    fit_unsorted <- fit_lgcp(locs = locs, sf = disc, smesh = smesh, tmesh = tmesh, parameters = sv)
    expect_equal(sum(fit_sorted$env$data$y), n)
    expect_equal(fit_sorted$env$data$y, fit_unsorted$env$data$y)
    ## a time knot with no events
    empty <- srt[srt$t < 1.5 | srt$t > 2.5, ]
    fit_empty <- fit_lgcp(locs = empty, sf = disc, smesh = smesh, tmesh = tmesh, parameters = sv)
    expect_equal(sum(fit_empty$env$data$y), nrow(empty))
    ## times outside the time mesh
    late <- srt
    late$t[nrow(late)] <- 4
    expect_error(fit_lgcp(locs = late, sf = disc, smesh = smesh, tmesh = tmesh, parameters = sv))
})

test_that("Spatiotemporal LGCP excludes mesh nodes outside the domain", {
    skip_on_cran()
    disc <- sf::st_sf(geometry = sf::st_geometry(sf::st_as_sf(spatstat.geom::disc())))
    smesh <- fmesher::fm_mesh_2d(loc.domain = sf::st_coordinates(disc)[, 1:2],
                                 max.edge = c(0.3, 0.6), offset = c(0.1, 0.5))
    tmesh <- fmesher::fm_mesh_1d(0:2)
    set.seed(2)
    r <- sqrt(runif(100)) * 0.9; th <- runif(100, 0, 2 * pi)
    locs <- data.frame(x = r * cos(th), y = r * sin(th), t = runif(100, 0, 2))
    fit <- fit_lgcp(locs = locs, sf = disc, smesh = smesh, tmesh = tmesh,
                    parameters = c(beta = 3, log_tau = 0, log_kappa = 0, atanh_rho = 0.2))
    expect_equal(fit$env$data$idx, as.numeric(fit$env$data$w > 0))
    expect_true(any(fit$env$data$idx == 0))
})

test_that("Simulate spatiotemporal LGCP with general time knots", {
    skip_on_cran()
    disc <- sf::st_sf(geometry = sf::st_geometry(sf::st_as_sf(spatstat.geom::disc())))
    smesh <- fmesher::fm_mesh_2d(boundary = fmesher::fm_as_segm(disc), max.edge = 0.3, cutoff = 0.1)
    for (knots in list(0:3, c(10, 12, 14))) {
        tmesh <- fmesher::fm_mesh_1d(knots)
        sim <- sim_lgcp(parameters = c(beta = 1, log_tau = 0, log_kappa = 0, atanh_rho = 0.5),
                        sf = disc, smesh = smesh, tmesh = tmesh)
        expect_equal(dim(sim$x), c(smesh$n, length(knots)))
        expect_true(all(is.finite(sim$y)))
    }
})

test_that("Marked LGCP: mark linked to the point process log intensity", {
    skip_on_cran()
    loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
    domain <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(loc.d))))
    smesh <- fmesher::fm_mesh_2d(loc.domain = loc.d, offset = c(0.3, 1),
                                 max.edge = c(0.3, 0.7), cutoff = 0.05)
    set.seed(3)
    sim <- sim_lgcp(parameters = c(beta = 4, log_tau = log(0.5), log_kappa = log(2)),
                    sf = domain, smesh = smesh)
    node <- rep(seq_len(smesh$n), times = sim$y)
    locs <- data.frame(x = smesh$loc[node, 1] + runif(length(node), -0.01, 0.01),
                       y = smesh$loc[node, 2] + runif(length(node), -0.01, 0.01))
    inside <- locs$x > 0 & locs$x < 3 & locs$y > 0 & locs$y < 3
    locs <- locs[inside, ]; node <- node[inside]
    ## mark = 5 + 0.5 * log intensity of the point process + noise (no mark field)
    m <- 5 + 0.5 * (4 + sim$x[node]) + rnorm(length(node), 0, 0.2)
    fit <- fit_mlgcp(locs = locs, marks = cbind(m = m), sf = domain, smesh = smesh,
                     ## start the mark intercept at the mean mark: starting at 0 lets nlminb stall on a
                     ## near-flat ridge (large marks_coefs_pp, near-zero field) on some platforms
                     parameters = list(betapp = 4, betamarks = matrix(mean(m), 1, 1), marks_coefs_pp = 0,
                                       log_tau = 0, log_kappa = 0),
                     methods = 0, fields = 0,
                     strfixed = matrix(log(0.2), nrow = nrow(locs), ncol = 1))
    cf <- get_coefs(fit)
    ## true values within 2.5 standard errors
    expect_lt(abs(cf["marks_coefs_pp", 1] - 0.5), 2.5 * cf["marks_coefs_pp", 2])
    expect_lt(abs(cf["betamarks", 1] - 5), 2.5 * cf["betamarks", 2])
})

test_that("fit_stelfi accepts mu and kernel volume is close to 1 on a fine mesh", {
    skip_on_cran()
    data(xyt, package = "stelfi")
    N <- 50
    locs <- data.frame(x = xyt$x[1:N], y = xyt$y[1:N])
    times <- xyt$t[1:N]
    domain <- sf::st_as_sf(xyt$window)
    bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    smesh <- fmesher::fm_mesh_2d(boundary = bnd, max.edge = 0.75, cutoff = 0.3)
    fit <- suppressWarnings(fit_stelfi(times = times, locs = locs, sf = domain, smesh = smesh,
                      parameters = list(mu = 3, alpha = 1, beta = 3, xsigma = 0.2,
                                        ysigma = 0.2, rho = 0.8)))
    expect_equal(fit$env$parameters$coefs, log(3))
    ## fine mesh over a square: kernel at the centre integrates to ~1
    sq <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(cbind(c(0, 4, 4, 0, 0), c(0, 0, 4, 4, 0))))))
    fine <- fmesher::fm_mesh_2d(boundary = fmesher::fm_as_segm(sq), max.edge = 0.05)
    w <- get_weights(fine, sq)$weights
    v <- stelfi:::kernel_mesh_volume(cbind(2, 2), fine$loc[, 1:2], w, 0.3, 0.2, 0.5)
    expect_equal(v, 1, tolerance = 0.02)
})
