#' delta hat squared for Chaiptak 2013 (helper function)
#'
#' @param ahatstar4 see function
#' @param p dimension
#' @param ahat2 see function
#' @param n sample size
#'
#' @keywords internal
#'
#'
deltahat2_func <- function(ahatstar4, p , ahat2, n){
  4 * ((2 * ahatstar4) / (p * (ahat2 ^ 2) * (n - 1)) + 1 / (n - 1) ^ 2)
}












