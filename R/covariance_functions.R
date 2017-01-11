#' Overall Covariance Matirx (helper function)
#'
#' @param A_ls sum of squares for group 1
#' @param ns sample size for group 1
#'
#' @keywords internal
#'
#'
overall_cov_func <- function(A_ls, ns){
  Reduce(`+`, A_ls) / Reduce(`+`, lapply(ns, function(x){x - 1}))
}

#' Sum of Squares (helper function)
#'
#' @param matrix data matrix
#'
#' @importFrom stats cov
#'
#' @keywords internal
#'
#'
A_func <- function(matrix){
  cov(matrix) * (nrow(matrix) - 1)
}

#' D diagonal matrix for Srivastava 2014 (helper funciton)
#'
#' @param matrix data matrix
#'
#' @keywords internal
#'
#'
Di_func <- function(matrix){
  diag(matrix %*% t(matrix))
}
