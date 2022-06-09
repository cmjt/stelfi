test_that("fit_hawkes estimation working", {
    compile_stelfi()
    # sim_hawkes example
    times <- sim_hawkes(mu = 0.3,alpha = 4, beta = 5)
    params <- c(mu = 0.3,alpha = 4, beta = 5)
    fit <- fit_hawkes(times = times,parameters = params)
    pars <- as.numeric(get_coefs(fit)[,1])
    expect_equal(pars[1], 0.2386207, tolerance = 1)
    expect_equal(pars[2], 1.8617510, tolerance = 1)
    expect_equal(pars[3], 3.5329454, tolerance = 1)
    
    # NIWA retweets example
    data(retweets_niwa)
    times <- unique(sort(as.numeric(difftime(retweets_niwa ,min(retweets_niwa),units = "mins"))))
    params <- c(mu = 9, alpha = 3, beta = 10)
    fit <- fit_hawkes(times = times, parameters = params) 

    pars <- as.numeric(get_coefs(fit)[,1])
    expect_equal(pars[1],0.06328, tolerance=1)
    expect_equal(pars[2],0.07597, tolerance=1)
    expect_equal(pars[3],0.07911, tolerance=1)
})
