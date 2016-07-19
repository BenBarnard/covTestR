#' Base for Brewer Kernel Quantile Estimator
#'
#' @param n sample size
#' @param i index
#' @param p quantile
#' @param X ith data point
#'
#' @keywords internal
#'
#' @export
#'
Bp_ <- function(n, i, p, X){
  ((n ^ (-1)) *
    ((2 * pi * p * (1 - p) / (n + 1)) ^ (-1 / 2)) *
    exp(- ((i / (n + 1) - p) ^ 2) / ((2 * p * (1 - p)) / (n + 1)))) * X
}

#' Brewer Kernel Quantile Estimator
#'
#' @param data vector of values
#' @param p quantile
#'
#' @keywords internal
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom plyr mlply
#'
#' @examples Bp_func(rnorm(100), .95)
Bp_func <- function(data, p){
  n <- length(data)
  i <- seq(1:n)
  X <- sort(data)
  Reduce(`+`, cbind(n, i, p, X) %>% mlply(Bp_))
}
