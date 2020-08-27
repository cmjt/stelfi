context("LGCP with INLA")

test_that("fit_lgcp_inla spatial", {
    win <- spatstat::owin(c(0, 3), c(0, 3))
    locs <- rlgcpspde(win, mu = 2, kappa = 1,sigma2 = 0.5)
    sp <- owin_to_sp(win)
    loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
    mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
                         max.edge = c(0.3, 0.7), cutoff = 0.05)
    fit <- fit_lgcp_inla(mesh = mesh,
                         locs = locs, sp = sp)
    expect_equal(fit$summary.hyperpar[2,1], 0.5, 0.2)
    expect_equal(fit$summary.fixed[1,1], 3, 1)
})

test_that("fit_lgcp_inla spatiotemporal", {
    win <- spatstat::owin(c(0, 3), c(0, 3))
    locs <- rlgcpspde(win, mu = 2, kappa = 1,sigma2 = 0.5, n = 3, rho = 0.8)
    sp <- owin_to_sp(win)
    loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
    mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
                         max.edge = c(0.3, 0.7), cutoff = 0.05)
    fit <- fit_lgcp_inla(mesh = mesh,
                         locs = locs[,1:2],temp = locs[,3], sp = sp)
    expect_equal(fit$summary.hyperpar[2,1], 0.5, 0.3)
    expect_equal(fit$summary.fixed[1,1], 3, 0.5)
})
