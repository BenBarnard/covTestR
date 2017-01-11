#' Tau for Chaipitak 2013 (helper function)
#'
#' @param ns sample size
#'
#' @keywords internal
#'
#'
tau_func <- function(ns){
  n <- Reduce(`+`, lapply(ns, function(x){x - 1}))
  ((n ^ 5) * ((n ^ 2) + n + 2)) /
    ((n + 1) * (n + 2) * (n + 4) * (n + 6) * (n - 1) * (n - 2) * (n - 3))
}
