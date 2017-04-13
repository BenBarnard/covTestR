#' Covariance Maker
#'
#' @param keepers keepers
#' @param offs offs
#' @param losers losers
#'
#' @importFrom stats toeplitz
#'
#' @export
#' @keywords internal
#' @examples cov_maker(keepers = list(c(20, 1, rep(5, 8))),
#'                     offs = list(.1, nrow = 10, ncol = 90),
#'                     losers = list(c(1, rep(.1, 89))))
cov_maker <- function(keepers, offs, losers){
  keep <- do.call(toeplitz, keepers)
  off <- do.call(matrix, offs)
  lose <- do.call(toeplitz, losers)
  rbind(cbind(keep, off), cbind(t(off), lose))
}
