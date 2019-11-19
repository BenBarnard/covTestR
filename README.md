
# covTestR

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/covTestR)](https://cran.r-project.org/package=covTestR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/covTestR)](https://cran.r-project.org/package=covTestR)
[![codecov](https://codecov.io/gh/BenBarnard/covTestR/branch/master/graph/badge.svg)](https://codecov.io/gh/BenBarnard/covTestR)
[![Travis-CI Build
Status](https://travis-ci.org/BenBarnard/covTestR.svg?branch=master)](https://travis-ci.org/BenBarnard/covTestR)

## Overview

covTestR is an equality of covariance testing suite. Currently only
equality of 2 and k group tests are available for high dimensional data.
There are future plans for one sample tests for high dimensinal data,
pretty print methods, and a wrapper function just to test.

In developing covTests we found it useful to have it play nice with the
“tidyverse.” At this point we have not thought of all the uses and
combinations with these packages so if you think of something not
currently implemented please file a minimal reproducible example on
github.

## Installation

You can install the latest development version from github with

``` r
if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
}
devtools::install_github("benbarnard/covTestR")
```

If you encounter a clear bug, please file a minimal reproducible example
on github.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tidyverse)
library(covTestR)

iris %>% group_by(Species) %>% homogeneityCovariances
```
