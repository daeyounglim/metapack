## Resubmission
This is a resubmission. In this version I have:

* removed unused functions in loglik_POCov.cpp and linearalgebra.cpp
* changed rwish to arma::wishrnd and deleted random.cpp
* created function `bmeta_analyze()` with a formula interface
* added `coef()` method
* renamed variables in the two data sets, and redocumented the data sets accordingly

## Test environments
* local x64 Windows 10 (Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz) R version 4.0.2 
* local macOS Big Sur 11.1 R 4.1.0
* ubuntu 20.04 (release and devel on R CMD check)
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel, release, and oldrelease)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.



## Downstream dependencies
There are no downstream dependencies for this package.

