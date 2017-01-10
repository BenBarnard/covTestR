#' Estimator for Ksi Srivastava and Yanagihara 2010 (helper function)
#'
#' @param n sample size
#' @param p dimension
#' @param ahat1 see function
#' @param ahat2 see function
#' @param ahat3 see function
#' @param ahat4 see function
#'
#' @keywords internal
#'
#'
ksihat2i_func <- function(n, p, ahat1, ahat2, ahat3, ahat4){
  (4 / ((n - 1) ^ 2)) * (((ahat2 ^ 2) / (ahat1 ^ 4)) +
                           ((2 * (n - 1)) / p) *
                           (((ahat2 ^ 3) / (ahat1 ^ 6)) -
                              ((2 * ahat2 * ahat3) / (ahat1 ^ 5)) +
                              (ahat4 / (ahat1 ^ 4))))
}
