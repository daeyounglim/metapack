## Resubmission
This is a resubmission. In this version I have:

* Removed the space in the doi specification to make it clickable
* Added executable examples for functions bayes.parobs and bayes.nmr
* OpenMP-related functions (model.comp.bayes.parobs and model.comp.bayesnmr) now default to two cores or fewer unless a user explicitly passes the # of cores to the functions
* Added two contributors that were missing to the package
* Description field in the DESCRIPTION file has been updated (minor changes)

## Test environments
* local macOS Big Sur 11.1 R 4.0.3
* ubuntu 20.04 (release and devel on R CMD check)
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, or WARNINGs.

There was 1 NOTEs (This is my first submission):

* checking CRAN incoming feasibility ... NOTE

The second note appears to be a false flag since the package does have a


## Downstream dependencies
There are currently no downstream dependencies for this package.

