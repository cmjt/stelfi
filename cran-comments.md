##  New submission

## Test environments

* ubuntu-latest (release & devel)
* local ubuntu 20.04
* windows-latest (release)
* macos-latest (release)
* macos-15-intel (release)

## R CMD check results

* No ERRORs or WARNINGs
 
The following NOTEs were received

New maintainer:
  Charlotte M. Jones-Todd <cmjt@st-andrews.ac.uk>
Old maintainer(s):
  Charlotte M. Jones-Todd <c.jonestodd@auckland.ac.nz>

 * I (Charlotte) am the maintainer and have recently moved institutions from the University of Auckland to the University of St Andrews and my email address has therefore changed accordingly

The following INFO was received

* checking installed package size ... INFO
  installed size is 91.5Mb
  sub-directories of 1Mb or more:
    libs  90.3Mb

 * This sub-directory is created when the C++ template, required for its functionality, is compiled.

## Downstream dependencies

* There are currently no downstream dependencies for this package.

## Changes

Maintainer email changed to cmjt@st-andrews.ac.uk

`fit_mhawkes()` (multivariate Hawkes model) updated coefficient names to include dimension when `get_coefs()` is called on fitted model object

Updated & added some thestthat tests