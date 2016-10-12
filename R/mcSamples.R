#' Simulate Monte Carlo Samples from a Multivariate Normal Distribution
#'
#' @param meanVec Numeric vector of means for the variables.
#' @param covMat Positive-definite symmetric matrix of the variables.
#' @param samples
#' @param pops
#' @param ... Other variables used in mvnorm function in the mass package.
#' @param matrix
#' @param tidy
#'
#' @return Matrix of maxN multivariate normal samples.
#'
#' @importFrom plyr mlply
#' @importFrom plyr ldply
#' @importFrom plyr llply
#' @importFrom MASS mvrnorm
#' @importFrom stats setNames
#' @importFrom reshape2 melt
#'
#'
#' @export
#'
#' @examples
#'
mcSamples <- function(meanVec, covMat, samples, pops, ..., matrix = FALSE, tidy = FALSE){
  pop_list <- pop_lists(meanVec, covMat, samples, pops)
  if(matrix == TRUE){
    llply(pop_list, function(list){
      mvrnorm(n = list$samples, mu = list$meanVec, Sigma = list$covMat)
    })
  }else{
    if(tidy == FALSE){
      ldply(pop_list, function(list){
        mvrnorm(n = list$samples, mu = list$meanVec, Sigma = list$covMat)
      }, .id = "population")
    }else{
      melt(ldply(pop_list, function(list){
        cbind(as.data.frame(mvrnorm(n = list$samples, mu = list$meanVec, Sigma = list$covMat)),
              samples = 1:samples)
      }, .id = "population"), id = c("population", "samples"))
    }
  }
}

#' List population parameters
#'
#' @param meanVec
#' @param covMat
#' @param samples
#' @param pops
#'
#' @return
#' @export
#'
#' @examples
pop_lists <- function(meanVec, covMat, samples, pops){
  params <- llply(seq(pops), function(pop, meanVec, covMat, samples){
    list(meanVec = if(is.list(meanVec)){meanVec[pop]}else{meanVec},
         covMat = if(is.list(covMat)){covMat[pop]}else{covMat},
         samples = if(is.list(samples)){samples[pop]}else{samples})
  }, meanVec = meanVec, covMat = covMat, samples = samples)
  setNames(params, 1:length(params))
}

