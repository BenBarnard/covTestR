setwd("~/Documents/R/Dissertation/critical_value")


library(EqualCov)

tests_ <- function(x, test, type){

  file <- readr::read_csv(paste0(type, "/reduction/", x))

  samples <- nrow(file) / (max(file$replication) * max(file$population))
  dimensions <- ncol(file[,-c(1,2)])

  readr::write_csv(plyr::ddply(file, plyr::.(replication),
                               function(file){
                                 do.call(test, list(x = file[,-1],
                                                    group = lazyeval::expr_find(population),
                                                    tidy = FALSE))
                               }),
                   paste(type, "/", "tests", "/", test, "_", x, sep = ""))
}

tests <- function(test, type){
  plyr::l_ply(list.files(paste0(type, "/reduction")), tests_, test = test, type = type)
}

tests("Chaipitak2013_test", "toeplitz")
tests("Schott2007_test", "toeplitz")
tests("Srivastava2007_test", "toeplitz")
tests("Srivastava2014_test", "toeplitz")
tests("SrivastavaYanagihara2010_test", "toeplitz")
tests("Chaipitak2013_test", "identity")
tests("Schott2007_test", "identity")
tests("Srivastava2007_test", "identity")
tests("Srivastava2014_test", "identity")
tests("SrivastavaYanagihara2010_test", "identity")
tests("diff_trace", "toeplitz")
tests("diff_trace", "identity")
tests("Ishii2016_test", "toeplitz")
tests("Ishii2016_test", "identity")
