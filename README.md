[![Build Status](https://travis-ci.org/daeyounglim/metapack.svg?branch=master)](https://travis-ci.org/daeyounglim/metapack)
[![R build status](https://github.com/daeyounglim/metapack/workflows/R-CMD-check/badge.svg)](https://github.com/daeyounglim/metapack/actions)
[![](https://img.shields.io/github/last-commit/daeyounglim/metapack.svg)](https://github.com/daeyounglim/metapack/commits/master)
[![](https://img.shields.io/github/languages/code-size/daeyounglim/metapack.svg)](https://github.com/daeyounglim/metapack)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


# metapack
**metapack** (version 0.1.0) provides functions to perform Bayesian inference for models introduced in the following papers:
+ Yao, H., Kim, S., Chen, M. H., Ibrahim, J. G., Shah, A. K., & Lin, J. (2015). Bayesian inference for multivariate meta-regression with a partially observed within-study sample covariance matrix. *Journal of the American Statistical Association*, **110(510)**, 528-544.
+ Li, H., Lim, D., Chen, M. H., Ibrahim, J. G., Kim, S., Shah, A. K., Lin, J. (2021). Bayesian network meta-regression hierarchical models using heavy-tailed multivariate random effects with covariate-dependent variances. Submitted (Under revision).

To see the model specification, please refer to the corresponding papers or the long-form vignette of this package.

## Installation
```r
# Currently, the package is available through GitHub
devtools::install_github("daeyounglim/metapack")
```

## Getting help
If you encounter a clear bug, please file an issue with a minimal reproducible example on [GitHub](https://github.com/daeyounglim/metapack/issues). For questions and other discussion, please email the [maintainer](mailto:daeyoung.lim@uconn.edu).

## Authors
+ Daeyoung Lim <daeyoung.lim@uconn.edu>
+ Ming-Hui Chen <ming-hui.chen@uconn.edu>
+ Sungduk Kim <kims2@mail.nih.gov>
+ Joseph G. Ibrahim <ibrahim@bios.unc.edu>

## Acknowledgments
+ Dr. Chen and Dr. Ibrahim's research was partially supported by NIH grants #GM70335 and #P01CA142538, and Merck & Co., Inc., Kenilworth, NJ, USA.
+ Dr. Kim's research was supported by the Intramural Research Program of National Institutes of Health, National Cancer Institute.