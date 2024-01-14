## Resubmission
This is a resubmission. In this version I have:

* Removed the C++11 requirement to accommodate updates in the BH package
* Updated the maintainer's e-mail address
* Added a citation file

## Test environments
* local macOS Sonoma 14.2.1 R 4.3.2
* Microsoft Windows Server 2019 (10.0.17763) on R CMD check
* Mac OS X (10.15.7) on R CMD check
* ubuntu 20.04 (release on R CMD check)
* win-builder (devel, release, and oldrelease)

## R CMD check results
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Daeyoung Lim <Daeyoung.Lim@fda.hhs.gov>'

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
