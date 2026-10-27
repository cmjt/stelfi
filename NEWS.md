# stelfi 1.0.4

## Breaking changes

 + `fit_mhawkes()`: the decay of the excitement of stream j is now governed by the decay rate of stream j (`beta[j]`), as documented, i.e., the intensity of stream j is mu_j + sum_k alpha_jk sum_{i in k} exp(-beta_j (t - t_i)). Previously the decay rate of the stream causing the excitement was used. `alpha[j, k]` is now constrained to be at most `beta[j]`. Estimates for models with unequal decay rates will change.
 + `fit_mlgcp()`: the coefficient linking each mark to the point process (`marks_coefs_pp`) now multiplies the log intensity of the point process. Previously it multiplied the log of the expected number of events in the area around each mesh node, which made the mark intercepts and links depend on the mesh. In addition, a mark with `fields = 0` no longer has the point process field added to its linear predictor. Estimates of `betamarks` and `marks_coefs_pp` will change.

## Minor improvements and fixes

 + `fit_lgcp()` (spatiotemporal): mesh nodes with no area within the domain are now excluded from the likelihood (as for the spatial model); previously they were included, biasing the estimates when the mesh extended beyond the domain.
 + `fit_lgcp()` (spatiotemporal): counts are now ordered by time knot regardless of the order of the data (previously unsorted data gave incorrect results), time knots with no events no longer cause a crash, events at the last knot are no longer dropped, and an informative error is given if any time is outside the range of the time knots.
 + `sim_lgcp()` (spatiotemporal): now works for any time mesh.
 + `fit_mhawkes()`: corrected the compensator, which subtracted the final event from every stream rather than only from the stream it belongs to.
 + `show_multivariate_hawkes()`: the fitted intensity now uses the correct (previously transposed) `alpha` matrix.
 + `fit_mlgcp()`: `locs` can now be a `data.frame` as well as a `matrix`.
 + `fit_stelfi()`: a starting value for the background rate can be given as `mu` (previously ignored unless given as `coefs`).
 + `fit_stelfi()`: warns if the mesh is too coarse to accurately integrate the estimated spatial self-excitation kernel over the domain.

# stelfi 1.0.3

## Breaking changes

 + None

## New features

 + None

## Minor improvements and fixes

 + Changed maintatiner email address to cmjt@st-andrews.ac.uk (previous c.jonestodd@auckland.ac.nz) due to moving institutions
 + For `fit_mhawkes()` (multivariate Hawkes model) have updated coefficient names to include dimension
 + Updated some tests to test for output structure

# stelfi 1.0.2

## Breaking changes

 + None

## New features

 + None

## Minor improvements and fixes

 + fixing CRAN NOTE for CRAN (re)release 


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
