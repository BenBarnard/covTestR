#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
#' @export
#'
ci_func <- function(n, p){
  p / n
}

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
#' @export
#'
c0_func <- function(n){
  n * ((n ^ 3) + 6 * (n ^ 2) + 21 * n + 18)
}

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
#' @export
#'
c1_func <- function(n){
  2 * n * (2 * (n ^ 2) + 6 * n + 9)
}

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
#' @export
#'
c2_func <- function(n){
  2 * n * (3 * n + 2)
}

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
#' @export
#'
c3_func <- function(n){
  n * (2 * (n ^ 2) + 5 * n + 7)
}
