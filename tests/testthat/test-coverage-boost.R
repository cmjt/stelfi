test_that("sim_hawkes rejects invalid input and method arguments", {
  expect_error(
    sim_hawkes(mu = 1, alpha = 2, beta = 1.5),
    "alpha must be less then beta"
  )

  expect_error(
    sim_hawkes(mu = 1, alpha = 0.5, beta = 2, method = "3"),
    "method should be either"
  )
})

test_that("fit_hawkes validates time, mark, and parameter constraints", {
  times <- c(0, 1, 2)

  expect_error(
    fit_hawkes(times = c(0, 0.5, 0.5), parameters = c(mu = 1, alpha = 0.3, beta = 0.8)),
    "simultaneous events"
  )

  expect_error(
    fit_hawkes(times = times, parameters = c(mu = 1, alpha = 0.9, beta = 0.8)),
    "alpha must be smaller"
  )

  expect_error(
    fit_hawkes(times = times, parameters = c(mu = 1, alpha = -0.2, beta = 0.8)),
    "alpha must be non-negative"
  )

  expect_error(
    fit_hawkes(times = times,
               parameters = c(mu = 1, alpha = 0.3, beta = 0.8),
               marks = c(1, -1, 1)),
    "marks cannot be negative"
  )
})

test_that("show_hawkes and show_hawkes_GOF exercise plotting and return branches", {
  data(retweets_niwa, package = "stelfi")
  times <- unique(sort(as.numeric(difftime(retweets_niwa,
                                          min(retweets_niwa),
                                          units = "mins"))))
  params <- c(mu = 9, alpha = 3, beta = 10)

  expect_error(
    show_hawkes(list(times = times, params = params), type = "bad"),
    "type.*fitted.*data.*both"
  )

  expect_s3_class(show_hawkes(list(times = times, params = params), type = "data"), "ggplot")
  expect_s3_class(show_hawkes(list(times = times, params = params), type = "both"), "gtable")

  gof <- show_hawkes_GOF(list(times = times, params = params),
                         plot = FALSE,
                         return_values = TRUE,
                         tests = FALSE)
  expect_is(gof, "list")
  expect_length(gof$compensator_differences, length(times) - 1)
})

test_that("get_weights and show_lambda hit mesh utility branches", {
  data(horse_mesh, package = "stelfi")
  data(horse_sf, package = "stelfi")

  weights <- get_weights(horse_mesh, horse_sf)
  expect_s3_class(weights, "sf")
  expect_true(all(weights$weights >= 0))

  if (requireNamespace("fmesher")) {
    data(xyt, package = "stelfi")
    domain <- sf::st_as_sf(xyt$window)
    locs <- data.frame(x = xyt$x, y = xyt$y)
    bnd <- fmesher::fm_as_segm(as.matrix(sf::st_coordinates(domain)[, 1:2]))
    smesh <- fmesher::fm_mesh_2d(boundary = bnd,
                                max.edge = 0.75,
                                cutoff = 0.3)

    fit <- fit_lgcp(locs = locs,
                    sf = domain,
                    smesh = smesh,
                    parameters = c(beta = 0, log_tau = log(1), log_kappa = log(1)))

    fields <- get_fields(fit, smesh)
    expect_true(is.matrix(fields) || is.vector(fields))
    expect_true(length(fields) > 0)

    plt <- show_lambda(fit, smesh = smesh, sf = domain)
    expect_s3_class(plt, "ggplot")
  }
})
