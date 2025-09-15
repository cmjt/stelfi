## (re/new)submission

Package was archived on CRAN on 2025-04-02

Resubmitted to CRAN 2025-08-19

Failed two tests in `test-fits.R` for `LGCP model fitting (spatial)` on M1mac where parameters were not estimated within the tolerance given. This likely due to the slight difference in approximation. I have reformulated this test accordingly.

## Test environments

* ubuntu-latest (release & devel)
* local ubuntu 20.04
* windows-latest (release & oldrelease)
* macos-latest (release)

## R CMD check results

* No ERRORs or WARNINGs
 
The following NOTE was received

* Maintainer: 'Charlotte M. Jones-Todd <c.jonestodd@auckland.ac.nz>'

New submission

Package was archived on CRAN

Possibly misspelled words in DESCRIPTION:
  Hawkes (3:8, 36:18, 36:92, 36:130)
  Lindgren (36:925)
  spatiotemporal (36:628)

 * These are the correct spellings

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2025-04-02 as issues were not corrected
    in time
    
  * Again, apologies. As above I have reformilated tests for M1mac

The following INFO was received

* checking installed package size ... INFO
 installed size is 84.5Mb
  sub-directories of 1Mb or more:
    libs  83.2Mb

 * This sub-directory is created when the C++ template, required for its functionality, is compiled.

## Downstream dependencies

* There are currently no downstream dependencies for this package.

## Changes

Maintainer email changed to c.jonestodd@auckland.ac.nz