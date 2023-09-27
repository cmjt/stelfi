# stelfi 0.0.3

## Breaking changes

 + None

## New features

 + Fitting of multivariate Hawkes processes `fit_mhawkes()`
 + Moved away from using `INLA` functionality to `fmesher`

## Minor improvements and fixes

 + Manually creating alias for package-level documentation
 + Additional testthat tests for fitting functions `fit_mhawkes()` and  `fit_mlgcp()`

# stelfi 0.0.2

## Breaking changes

 + None

## New features

 + None

## Minor improvements and fixes

 + Bug fix for `show_hawkes_GOF()`, compensator plotted on goodness of fit plots
 + Additional testthat tests for utility functions `meshmetrics()`, `show_hawkes_GOF()`, and `get_fields()`
 + Additional testthat tests for fitting/simulation functions `sim_hawkes(..., method = "2")` and `sim_lgcp()`


# stelfi 0.0.1

## Initial release

Model fitting functions include

 + `fit_hawkes()` and `fit_hawkes_cbf()` to fit self-exciting Hawkes processes to temporal point pattern data
 + `fit_lgcp()` fits a log-Gaussian Cox process to either spatial or spatiotemporal point pattern data. If a spatiotemporal model is fitted an AR1 process is assumed for the temporal progression
 + `fit_mlgcp()` fits a joint likelihood model between the point locations and mark(s)
 + `fit_stelfi()` fits self-exciting spatiotemporal Hawkes models to point pattern data. The self-excitement is Gaussian in space and exponentially decaying over time. In addition, a GMRF can be included to account for latent spatial dependency
