## Resubmission
This is a resubmission. In this version I have:

* Wishart samplers have been coded independently because RcppArmadillo's seeding does not guarantee reproducibility
* dots in non-S3 functions have been replaced with underscores, which is a minor change since these functions are not public-facing
* changed the DOI link to the journal webpage link in a function documentation to avoid the DOI check

## Test environments
* local macOS Big Sur 11.1 R 4.1.3
* Microsoft Windows Server 2019 (10.0.17763) on R CMD check
* Mac OS X (10.15.7) on R CMD check
* ubuntu 20.04 (release on R CMD check)
* win-builder (devel, release, and oldrelease)

## R CMD check results
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Daeyoung Lim <daeyoung.lim@uconn.edu>'

Found the following (possibly) invalid URLs:
  URL: https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1044
    From: man/bmeta_analyze.Rd
    Status: 503
    Message: Service Unavailable
Found the following (possibly) invalid DOIs:
  DOI: 10.1002/sim.8983
    From: DESCRIPTION
    Status: Service Unavailable

    Message: 503

The URL and DOI are valid and I've checked that they work, but I couldn't figure out a way to make this note go away.

## Downstream dependencies
There are no downstream dependencies for this package.
