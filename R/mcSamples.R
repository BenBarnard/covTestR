#' Simulate Monte Carlo Samples from a Multivariate Normal Distribution
#'
#' @param meanVec Numeric vector of means for the variables.
#' @param covMat Positive-definite symmetric matrix of the variables.
#' @param maxN Numeric scalar of the maximum of the sample sizes in the simulation.
#' @param ... Other variables used in mvnorm function in the mass package.
#'
#' @return Matrix of maxN multivariate normal samples.
#'
#' @importFrom dplyr rename
#' @importFrom plyr ldply
#' @importFrom stats setNames
#' @importFrom reshape2 melt
#' @importFrom MASS mvrnorm
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' mcSamples(c(0,0,0), diag(1, 3), 10, 2)
mcSamples <- function(meanVec, covMat, maxN, pops, ...){
  replicate(pops,
            mvrnorm(n = maxN, mu = meanVec, Sigma = covMat) %>%
              melt %>%
              setNames(c('Subjects', 'Variables', 'Value')),
            simplify = FALSE) %>%
    setNames(c(1:pops)) %>%
    ldply %>%
    rename(Group = `.id`)
}
