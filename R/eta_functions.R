#' Estimator for eta squared Srivastave 2007 (helper function)
#'
#' @param n sample size
#' @param p dimension
#' @param ahat4 see function
#' @param ahat2 see function
#'
#' @keywords internal
#'
#' @export
#'
etahat2i_func <- function(n, p, ahat4, ahat2){
  (4 / ((n - 1) ^ 2)) *
    (ahat2 ^ 2) *
    (1 + ((2 * (n - 1) * ahat4) / (p * ahat2 ^ 2)))
}
