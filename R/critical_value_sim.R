#' Sim Function to get data for Critical Value Data
#'
#' @param data the data
#' @param ... other stuff
#'
#'
#' @export
#'
#' @keywords internal
#'
#' @examples critical_value_sim(list(mean = c(0,0,0),
#'                                   cov = diag(1, 3),
#'                                   samples = 10,
#'                                   populations = 2),
#'                              10000, .95, Chaipitak2013_test)
critical_value_sim <- function(params, replications, quantile, test){

  replicate(replications)
}


