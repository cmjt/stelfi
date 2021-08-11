test_that("fit_hawkes estimation working", {
    times <- sim_hawkes(mu = 0.3,alpha = 4, beta = 5)
    params <- c(mu = 0.3,alpha = 4, beta = 5)
    fit <- fit_hawkes(times = times,parameters = params)
    pars <- as.numeric(get_coefs(fit)[,1])
    expect_equal(pars[1], 0.2386207, tolerance = 1)
    expect_equal(pars[2], 1.8617510, tolerance = 1)
    expect_equal(pars[3], 3.5329454, tolerance = 1)
})
