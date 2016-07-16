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
  ((n ^ (-1)) * (gamma(n + 1) / ((gamma(i)) * gamma(n - i + 1))) * (p ^ (i - 1)) * ((1 - p) ^ (n - i))) * X
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
Bp_func <- function(data, p){
  browser()
  n <- length(data)
  data %<>% as.data.frame

}
