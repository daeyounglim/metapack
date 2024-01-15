## Resubmission
This is a resubmission. In this version I have:

* Removed the C++11 requirement to accommodate updates in the BH package
* Updated the maintainer's e-mail address
* Added a citation file

## Test environments
* local macOS Sonoma 14.2.1 R 4.3.2 R CMD check
* win-builder (devel, release, and oldrelease)
* R-hub builder (https://builder.r-hub.io)
    * Windows Server 2022, R-devel, 64 bit
        > ❯ checking CRAN incoming feasibility ... NOTE
            Maintainer: 'Daeyoung Lim <Daeyoung.Lim@fda.hhs.gov>'****
  
            New maintainer:
              Daeyoung Lim <Daeyoung.Lim@fda.hhs.gov>
            Old maintainer(s):
              Daeyoung Lim <daeyoung.lim@uconn.edu>

        > ❯ checking HTML version of manual ... NOTE
            Skipping checking math rendering: package 'V8' unavailable

        > ❯ checking for non-standard things in the check directory ... NOTE
            Found the following files/directories:
              ''NULL''

        > ❯ checking for detritus in the temp directory ... NOTE
            Found the following files/directories:
            'lastMiKTeXException'
    * Ubuntu Linux 20.04.1 LTS, R-release, GCC
    * Fedora Linux, R-devel, clang, gfortran
    * Debian Linux, R-devel, GCC ASAN/UBSAN

## Downstream dependencies
There are no downstream dependencies for this package.
