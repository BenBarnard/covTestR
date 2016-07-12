#' Tau for Chaipitak 2013 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
#'
tau_func <- function(n1, n2){
  n <- n1 + n2 - 2
  ((n ^ 5) * ((n ^ 2) + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6) * (n - 1) * (n - 3))
}
