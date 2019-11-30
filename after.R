install.packages(c("pkgdown", "covr"))

library(pkgdown)
library(covr)

pkgdown::build_site()
covr::codecov()