## Resubmission

This is a resubmission. In this version I have fixed the following issues as requested

* \value tags added to all required exported functions
* Examples unwrapped if they execute in < 5s
* \dontrun replaced with \donttest
* Examples using Suggests packages wrapped in `requireNamespace`
* `old = options()` ; `on.exit(options(old))` added to `.onLoad` in `zzz.R` so that user settings are reset
when function is exited
* Ensures that no packages are installed within any function. `stelfi_load_inla` in `zzz.R` prompts user to do so by printing the commend required.

## Test environments

* Local Ubuntu 18.04 & 20.04 install (release)
* win-builder (devel and release)
* macos-builder (release)

## R CMD check results

* No ERRORs or WARNINGs

On all tests the following NOTEs were received

*  checking CRAN incoming feasibility ... NOTE
   Maintainer: 'Charlotte M. Jones-Todd <cmjonestodd@gmail.com>'

   New submission

   * This is the first package I have submitted to CRAN.

* checking installed package size ... NOTE
  installed size is 69.7Mb
  sub-directories of 1Mb or more:
    libs  68.6Mb

   * This sub-directory is created when the C++ template, required for its functionality, is compiled.

* checking package dependencies ... NOTE
Package suggested but not available for checking: ‘INLA’
* Suggests or Enhances not in mainstream repositories:
  INLA

   * INLA is a non-CRAN package, it is listed under `Additional_repositories`

On win-builder the additional an NOTE is received

* Possibly misspelled words in DESCRIPTION:
  Hawkes (3:8, 32:18, 32:92, 32:130)
  Lindgren (32:925)
  spatiotemporal (32:628)
  Found the following (possibly) invalid DOIs:
  DOI: 10.1111/j.1467-9868.2011.00777.x
    From: DESCRIPTION

    * spatiotemporal is the correct spelling; Hawkes is the correct spelling, named after the statistician Alan G. Hawkes; Lindgren is again the correct spelling, the surname of the researcher Finn Lindgren. I have manually checked the DOI and it is correct.
