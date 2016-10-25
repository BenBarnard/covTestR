setwd("~/Documents/R/Dissertation/power")

library(EqualCov)

crit_data <- function(type, dimensions, samples, difference, populations, replications){
  plyr::m_ply(expand.grid(dimensions = dimensions,
                          samples = samples,
                          difference = difference),
              function(dimensions, samples, difference, populations, replications, type){
                if(type == "identity"){
                  covMat <- list(diag(1, dimensions), diag(1 + sqrt(difference / dimensions), dimensions))
                }

                if(type == "toeplitz"){
                  mat <- toeplitz(.5 ^ seq(0, (dimensions - 1)))
                  covMat <- list(mat, mat + sqrt(difference / (dimensions ^ 2)))
                }

                readr::write_csv(plyr::rdply(.n = replications,
                                       mcSamples(meanVec = rep(0, nrow(covMat[[1]])),
                                                 covMat = covMat,
                                                 samples = samples,
                                                 pops = populations),
                                       .id = "replication"),
                                 paste(type, "/", difference, "_", dimensions, "_", samples, ".csv", sep = ""))
                },
              replications = replications,
              type = type, populations = populations)
}

crit_data(type = "identity",
          dimensions = c(20, 25, 30, 35, 40, 45, 50),
          samples = c(15),
          difference = c(1, 2, 3, 4, 5),
          populations = 2,
          replications = 10000)

crit_data(type = "toeplitz",
          dimensions = c(20, 25, 30, 35, 40, 45, 50),
          samples = c(15),
          difference = c(1, 2, 3, 4, 5),
          populations = 2,
          replications = 10000)



