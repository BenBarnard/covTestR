#' Sim Function to get data for Critical Value Data
#'
#' @param data the data
#' @param ... other stuff
#'
#' @importFrom magrittr %>%
#' @importFrom plyr laply
#'
#' @export
#'
#' @keywords internal
#'
#' @examples critical_value_sim(list(meanVec = c(0,0,0),
#'                                   covMat = diag(1, 3),
#'                                   samples = 10,
#'                                   pops = 2),
#'                              100, .95, Chaipitak2013_test)
critical_value_sim <- function(params, replications, quantile, test){
  replicate(replications, do.call(mcSamples, params), simplify = FALSE) %>%
    laply(test) %>%
    Bp_func(quantile)
}


