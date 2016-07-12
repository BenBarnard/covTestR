#' delta hat squared for Chaiptak 2013 (helper function)
#'
#' @param ahatstar4 see function
#' @param p dimension
#' @param ahat2 see function
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#'
#' @keywords internal
#'
deltahat2_func <- function(ahatstar4, p , ahat2, n1, n2){
  4 * (((2 * ahatstar4) / (p * ahat2 ^ 2)) * sum(1 / (n1 - 1), 1 / (n2 - 1)) + sum(1 / ((n1 - 1) ^ 2), 1 / ((n2 - 1) ^ 2)))
}












