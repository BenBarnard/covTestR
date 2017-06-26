#' Trace of Matrix
#'
#'
#'
#' @export
#' @keywords internal
tr <- function(x){
  sum(diag(x))
}
