setwd("~/Documents/R/Dissertation/critical_value")


library(EqualCov)
library(classReduce)

dim_reduce_ <- function(x, reduction, type, reducedDim){

  file <- readr::read_csv(paste0(type, "/", x))

  samples <- nrow(file) / (max(file$replication) * max(file$population))
  dimensions <- ncol(file[,-c(1,2)])

  readr::write_csv(plyr::ddply(file, plyr::.(replication),
                               function(file){
                                 do.call(reduction, list(x = file[,-1],
                                                    group = lazyeval::expr_find(population),
                                                    targetDim = reducedDim))$reducedData
                               }),
                   paste(type, "/", "reduction", "/", reduction, "_", dimensions, "_", samples, "_", reducedDim, ".csv", sep = ""))
}

reductions <- function(reduction, type, reducedDim){
  plyr::l_ply(list.files(paste0(type)), dim_reduce_, reduction = reduction, type = type, reducedDim = reducedDim)
}

reductions("SConcat", "identity", 20)
reductions("DataConcatScatter", "identity", 20)
reductions("DataConcatScatterBlock", "identity", 20)
reductions("SConcat", "identity", 19)
reductions("DataConcatScatter", "identity", 19)
reductions("DataConcatScatterBlock", "identity", 19)
reductions("SConcat", "identity", 18)
reductions("DataConcatScatter", "identity", 18)
reductions("DataConcatScatterBlock", "identity", 18)
reductions("SConcat", "identity", 17)
reductions("DataConcatScatter", "identity", 17)
reductions("DataConcatScatterBlock", "identity", 17)
reductions("SConcat", "identity", 16)
reductions("DataConcatScatter", "identity", 16)
reductions("DataConcatScatterBlock", "identity", 16)
reductions("SConcat", "identity", 15)
reductions("DataConcatScatter", "identity", 15)
reductions("DataConcatScatterBlock", "identity", 15)
reductions("SConcat", "identity", 14)
reductions("DataConcatScatter", "identity", 14)
reductions("DataConcatScatterBlock", "identity", 14)
reductions("SConcat", "toeplitz", 20)
reductions("DataConcatScatter", "toeplitz", 20)
reductions("DataConcatScatterBlock", "toeplitz", 20)
reductions("SConcat", "toeplitz", 19)
reductions("DataConcatScatter", "toeplitz", 19)
reductions("DataConcatScatterBlock", "toeplitz", 19)
reductions("SConcat", "toeplitz", 18)
reductions("DataConcatScatter", "toeplitz", 18)
reductions("DataConcatScatterBlock", "toeplitz", 18)
reductions("SConcat", "toeplitz", 17)
reductions("DataConcatScatter", "toeplitz", 17)
reductions("DataConcatScatterBlock", "toeplitz", 17)
reductions("SConcat", "toeplitz", 16)
reductions("DataConcatScatter", "toeplitz", 16)
reductions("DataConcatScatterBlock", "toeplitz", 16)
reductions("SConcat", "toeplitz", 15)
reductions("DataConcatScatter", "toeplitz", 15)
reductions("DataConcatScatterBlock", "toeplitz", 15)
reductions("SConcat", "toeplitz", 14)
reductions("DataConcatScatter", "toeplitz", 14)
reductions("DataConcatScatterBlock", "toeplitz", 14)
