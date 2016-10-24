setwd("~/Documents/R/Dissertation/critical_value")


crit_data <- function(type, dimensions, samples, populations, replications){
  plyr::m_ply(expand.grid(dimensions = dimensions,
                          samples = samples),
              function(dimensions, samples, populations, replications, type){

                if(type == "identity"){
                  covMat <- diag(1, dimensions)
                }

                if(type == "toeplitz"){
                  covMat <- toeplitz(.5 ^ seq(0, (dimensions - 1)))
                }

                readr::write_csv(plyr::rdply(.n = replications,
                                       mcSamples(meanVec = rep(0, nrow(covMat)),
                                                 covMat = covMat,
                                                 samples = samples,
                                                 pops = populations,
                                                 matrix = FALSE,
                                                 tidy = FALSE),
                                       .id = "replication"),
                                 paste(type, "/", dimensions, "_", samples, ".csv", sep = ""))
                },
              replications = replications,
              type = type, populations = populations)
}

crit_data(type = "",
          dimensions = c(20, 25, 30, 35, 40, 45, 50),
          samples = c(15),
          populations = 2,
          replications = 10000)

