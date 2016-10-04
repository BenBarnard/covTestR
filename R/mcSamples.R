#' Simulate Monte Carlo Samples from a Multivariate Normal Distribution
#'
#' @param meanVec Numeric vector of means for the variables.
#' @param covMat Positive-definite symmetric matrix of the variables.
#' @param maxN Numeric scalar of the maximum of the sample sizes in the simulation.
#' @param ... Other variables used in mvnorm function in the mass package.
#'
#' @return Matrix of maxN multivariate normal samples.
#'
#' @importFrom plyr mlply
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#'
mcSamples <- function(meanVec, covMat, samples, pops, ...){
  pop_list <- pop_lists(meanVec, covMat, samples, pops)
  llply(pop_list, function(list){
        mvrnorm(n = list$samples, mu = list$meanVec, Sigma = list$covMat)
    })
}

pop_lists <- function(meanVec, covMat, samples, pops){
  llply(seq(pops), function(pop, meanVec, covMat, samples){
    list(meanVec = if(is.list(meanVec)){meanVec[pop]}else{meanVec},
         covMat = if(is.list(covMat)){covMat[pop]}else{covMat},
         samples = if(is.list(samples)){samples[pop]}else{samples})
  }, meanVec = meanVec, covMat = covMat, samples = samples)
}

