#' Covariance Maker
#'
#' @param keepers
#' @param offs
#' @param losers
#'
#' @return
#' @export
#'
#' @examples cov_maker(keepers = list(c(20, 1, rep(5, 8))),
#'                     offs = list(0, nrow = 10, ncol = 90),
#'                     losers = list(c(1, rep(0, 89))))
cov_maker <- function(keepers, offs, losers){
  keep <- do.call(toeplitz, keepers)
  off <- do.call(matrix, offs)
  lose <- do.call(toeplitz, losers)
  rbind(cbind(keep, off), cbind(t(off), lose))
}
