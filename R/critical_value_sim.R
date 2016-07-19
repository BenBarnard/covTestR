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
#' @examples critical_value_sim_(mcSamples(c(0,0,0), diag(1, 3), 10, 2))
critical_value_sim_ <- function(data, quantile, ..., dots){

}

#' Critical Value Simulation
#'
#' @param sim_size
#' @param sim_data_func
#'
#' @export
#'
#' @importFrom lazyeval lazy
#' @importFrom lazyeval lazy_eval
#' @importFrom lazyeval lazy_dots
#' @importFrom magrittr %>%
#' @importFrom plyr llply
#'
#' @keywords internal
#'
critical_value_sim <- function(sim_size, sim_data_func, quantile, ...){
  lazydata <- lazy(sim_data_func)
  dots <- lazy_dots(...)
  replicate(sim_size, lazy_eval(lazydata), simplify = FALSE) %>% llply(critical_value_sim_, dots = dots, quantile = quantile)
}

