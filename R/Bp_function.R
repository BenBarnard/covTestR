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
     (gamma(n + 1) /
        ((gamma(i)) * gamma(n - i + 1))) *
     (p ^ (i - 1)) * ((1 - p) ^ (n - i))) * X
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
Bp_func <- function(data, p){
  n <- length(data)
  i <- seq(1:n)
  X <- sort(data)
  Reduce(`+`, cbind(n, i, p, X) %>% mlply(Bp_))
}
