## Test environments
* local macOS Big Sur 11.1 R 4.0.3
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

# macOS Big Sur 11.1
## 2.3 GHz Quad-Core Intel Core i7

── R CMD check results ─────────────────────────────────────── metapack 0.1 ────
Duration: 3m 58.7s

> checking R code for possible problems ... NOTE
  plot.sucra: no visible binding for global variable ‘trt’
  plot.sucra: no visible binding for global variable ‘CDF’
  plot.sucra: no visible binding for global variable ‘PDF’
  Undefined global functions or variables:
    CDF PDF trt

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded



# WinBuilder

* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking serialization versions ... OK
* checking whether package 'metapack' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking for future file timestamps ... OK
* checking 'build' directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* loading checks for arch 'i386'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* loading checks for arch 'x64'
** checking whether the package can be loaded ... OK
** checking whether the package can be loaded with stated dependencies ... OK
** checking whether the package can be unloaded cleanly ... OK
** checking whether the namespace can be loaded with stated dependencies ... OK
** checking whether the namespace can be unloaded cleanly ... OK
** checking loading without being on the library search path ... OK
** checking use of S3 registration ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... [14s] NOTE
plot.sucra: no visible binding for global variable 'trt'
plot.sucra: no visible binding for global variable 'CDF'
plot.sucra: no visible binding for global variable 'PDF'
Undefined global functions or variables:
  CDF PDF trt
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd line widths ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of 'data' directory ... OK
* checking data for non-ASCII characters ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking line endings in C/C++/Fortran sources/headers ... OK
* checking line endings in Makefiles ... OK
* checking compilation flags in Makevars ... OK
* checking for GNU extensions in Makefiles ... OK
* checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS) ... OK
* checking use of PKG_*FLAGS in Makefiles ... OK
* checking use of SHLIB_OPENMP_*FLAGS in Makefiles ... OK
* checking pragmas in C/C++ headers and code ... OK
* checking compiled code ... OK
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking examples ...
** running examples for arch 'i386' ... [2s] OK
** running examples for arch 'x64' ... [2s] OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes in 'inst/doc' ... OK
* checking re-building of vignette outputs ... [2s] OK
* checking PDF version of manual ... OK
* checking for detritus in the temp directory ... OK
* DONE
Status: 2 NOTEs