## Resubmission

This is a resubmission. In this version I have:

* Added a stelfi-package \alias, as requested by email from Kurt Hornik.
* Removed all dependency in the non-CRAN package INLA; instead functions from the new CRAN package fmesher are called.
* Added an additional model fitting function for a multivariate Hawkes process, as well as associated testthat tests.


## Test environments

* ubuntu-latest (release & devel)
* local ubuntu 18.04 & 20.04
* windows-latest (release)
* macos-latest (release)

## R CMD check results

* No ERRORs or WARNINGs

On all tests the following NOTE was received

* checking installed package size ... NOTE
 installed size is 62.7Mb 
 sub-directories of 1Mb or more:
 libs  61.5Mb

 * This sub-directory is created when the C++ template, required for its functionality, is compiled.

## Downstream dependencies

* There are currently no downstream dependencies for this package.