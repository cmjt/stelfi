context("self-exciting Hawkes process")

test_that("fit_hawkes works", {
    t <- sim_hawkes(mu = 3, alpha = 0.5, beta = 2, n = 100, plot = FALSE,
                    seed = 1222, method = "1")
    params <- c(mu = 3,alpha = 0.5,beta = 2)
    fit <- fit_hawkes(times = t,parameters = params,upper = c(50,5,15),lower = c(0,0,0))
    expect_equal(fit$par.fixed[1], 3, 0.5)
})
