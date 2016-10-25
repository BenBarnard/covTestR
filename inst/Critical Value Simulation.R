library(plyr)
library(tidyverse)
library(lazyeval)

critical_simulation <- function(data){

  tests <- c("Chaipitak2013_test", "Ishii2016_test", "Schott2007_test",
             "Srivastava2007_test", "Srivastava2014_test",
             "SrivastavaYanagihara2010", "boxM")

  svdMethod <- c("svd", "sparsesvd")

  reductions <- c("SConcat", "DataConcatScatter",
                  "DataConcatScatterBlock")


exgrid <- expand.grid(tests, reductions, svdMethod)
original <- expand.grid(tests[!(tests == "boxM")], "Original", "none")
situations <- rbind(exgrid, original)


}


critical_simulation(mcSamples(rep(0, 100), diag(1, 100), 15, 2))
