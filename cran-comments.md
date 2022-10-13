## Test environments

* Local Xubuntu 18.04 & 20.04 install (release)
* win-builder (devel and release)
* macos-builder (release)

## R CMD check results

* No ERRORs or WARNINGs

On all tests the following NOTEs were recieved

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
