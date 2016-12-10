#' Title
#'
#' @param type
#' @param dimensions
#' @param samples
#' @param difference
#' @param populations
#' @param replications
#' @param directory
#'
#' @return
#' @export
#'
#' @importFrom plyr m_ply
#' @importFrom plyr rdply
#' @importFrom plyr ddply
#' @importFrom plyr ldply
#' @importFrom lazyeval expr_find
#' @importFrom classReduce DataConcatScatter
#' @importFrom classReduce DataConcatScatterBlock
#' @importFrom classReduce SConcat
#' @importFrom classReduce SDiff
#' @importFrom readr write_csv
#'
#'
#' @examples crit_data3(type = "identity",
#'                     dimensions = c(100, 125, 150, 175, 200),
#'                     samples = c(15),
#'                     difference = c(0, 1, 2, 3, 4, 5),
#'                     populations = 3,
#'                     replications = 1000,
#'                     reductionMethods = c("DataConcatScatter", "DataConcatScatterBlock",
#'                     "SConcat", "SDiff"),
#'                     directory = "~/Documents/R/Dissertation/data")
#'
crit_data3 <- function(type, dimensions, samples, difference, populations, replications, reductionMethods, directory){
  comn <-  expand.grid(dimensions = dimensions,
                       samples = samples,
                       difference = difference)

  m_ply(comn,
        function(dimensions, samples, difference, populations, replications, type, reductionMethods, directory){
          if(type == "identity"){
            covMat <- list(diag(1, dimensions), diag(1, dimensions), diag(1 + difference, dimensions))
          }

          if(type == "toeplitz"){
            mat <- toeplitz(.5 ^ seq(0, (dimensions - 1)))
            covMat <- list(mat, mat, mat + difference)
          }

          if(type == "elliptical"){
            mat <- cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                             offs = list(0, nrow = 10, ncol = dimensions - 10),
                             losers = list(c(1, rep(0, dimensions - 11))))
            covMat <- list(mat, mat, mat + difference)
          }

          if(type == "elliptical2"){
            mat <- cov_maker(keepers = list(c(20, 1, rep(5, 8))),
                             offs = list(0, nrow = 10, ncol = dimensions - 10),
                             losers = list(c(1, rep(0, dimensions - 11))))
            covMat <- list(mat, mat * difference)
          }

          originaldata <- rdply(.n = replications,
                                mcSamples(meanVec = rep(0, nrow(covMat[[1]])),
                                          covMat = covMat,
                                          samples = samples,
                                          pops = populations),
                                .id = "replication")
          names(originaldata) <- c(names(originaldata)[names(originaldata) %in%
                                                         c("replication", "population")],
                                   paste0("V", names(originaldata)[!(names(originaldata) %in%
                                                                       c("replication", "population"))]))
          originaldata <- cbind(originaldata, data.frame(originaldimensions = dimensions,
                                                         difference = difference,
                                                         ReductionMethod = "None",
                                                         type = type))

          write_csv(rbind(ldply(lapply(reductionMethods, function(reduction){
            cbind(ldply(lapply(1:replications, function(replication){
              cbind(do.call(reduction, list(x = originaldata[originaldata$replication == replication,
                                                             !(names(originaldata) %in% c("originaldimensions",
                                                                                          "difference",
                                                                                          "replication",
                                                                                          "ReductionMethod",
                                                                                          "type"))],
                                            group = expr_find(population),
                                            targetDim = dimensions))$reducedData,
                    originaldata[originaldata$replication == replication,
                                 names(originaldata) %in% c("originaldimensions",
                                                            "difference",
                                                            "replication")]
              )
            })),
            data.frame(ReductionMethod = rep(reduction, replications * populations * samples),
                       type = rep(type, replications * populations * samples)))
          })), originaldata), paste0(directory, "/", dimensions, " ", samples,
                                     " ", difference, " ", type, " 3 .csv"))
        },
        replications = replications,
        type = type,
        populations = populations,
        reductionMethods = reductionMethods,
        directory = directory)
}
