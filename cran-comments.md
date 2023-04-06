## Resubmission

This is a resubmission. In this version I have:

* Corrected the Date field in the DESCRIPTION.
* Removed use of the previously undefined global functions density & dexp, which
were called internally within ggplot2 plotting code.


## Test environments

* ubuntu-latest (release & devel & oldrel-1)
* local ubuntu 18.04 & 20.04
* windows-latest (release)
* macos-latest (release)

## R CMD check results

* No ERRORs or WARNINGs

On all tests the following NOTE was received

Package suggested but not available for checking: ‘INLA’

   * INLA is a non-CRAN package, it is listed under `Additional_repositories`

