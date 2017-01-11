#' Estimator for Gamma Srivastava and Yanagihara 2010 (helper function)
#'
#' @param ahat2i see function
#' @param ahat1i see function
#'
#' @keywords internal
#'
#'
gammahati_func <- function(ahat2i, ahat1i){
  ahat2i / (ahat1i ^ 2)
}
