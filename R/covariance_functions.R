#' Overall Covariance Matirx (helper function)
#'
#' @param A1 sum of squares for group 1
#' @param A2 sum of squares for group 2
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#'
#' @keywords internal
#'
#' @export
#'
overall_cov_func <- function(A1, A2, n1, n2){
  (1 / (n1 + n2 - 2)) * (A1 + A2)
}

#' Sum of Squares (helper function)
#'
#' @param matrix data matrix
#'
#' @keywords internal
#'
#' @export
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
#' @export
#'
Di_func <- function(matrix){
  diag(matrix %*% t(matrix))
}
