#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param n sample size
#' @param p dimension
#' @param sample_cov sample covariance matrix
#'
#' @keywords internal
#'
ahat2i_func <- function(n, p, sample_covs){
  (((n - 1) ^ 2) / (p * (n -2) * (n + 1))) *
    (tr(sample_covs ^ 2) - (1 / (n - 1)) * (tr(sample_covs)) ^ 2)
}

#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param n sample size
#' @param p dimension
#' @param sample_cov sample covariance matrix
#'
#' @keywords internal
#'
ahat2_func <- function(n1, n2, p, overall_cov){
  n <- n1 + n2 - 2
  (((n) ^ 2) / (p * (n - 1) * (n + 2))) *
    (tr(overall_cov ^ 2) - (1 / (n)) * (tr(overall_cov)) ^ 2)
}

#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param p dimension
#' @param sample_cov sample covariance matrix
#'
#' @keywords internal
#'
ahat1_func <- ahat1i_func <- function(p, sample_cov){
  (1 / p) * tr(sample_cov)
}


#' Estimator for expansion term in Frobenius Norm Chaipitak 2013 (helper funciton)
#'
#' @param tau see function
#' @param p dimension
#' @param sample_cov sample covariance matrix
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#'
#' @aliases ahat1i_func
#'
#' @keywords internal
#'
ahatStar4_func <- function(tau, p, sample_cov, n1, n2){
  n <- n1 + n2 - 2
  (tau / p) *
    (tr(sample_cov ^ 4) +
       (-4 / n) * tr(sample_cov ^ 3) * tr(sample_cov) +
       (-((2 * n ^ 2) + 3 * n - 6) / (n * ((n ^ 2) + n + 2))) * (tr(sample_cov ^ 2) ^ 2) +
       ((2 * (5 * n + 6)) / (n * ((n ^ 2) + n + 2))) * tr(sample_cov ^ 2) * (tr(sample_cov) ^ 2) +
       (-(5 * n + 6) / ((n ^ 2) * ((n ^ 2) + n + 2))) * tr(sample_cov) ^ 4)
}

#' Estimator for frobenius norm expansion term Srivastava 2007
#'
#' @param A1 sum of squares for group 1
#' @param A2 sum of squares for group 2
#' @param p dimension
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#' @param ahat2 see function
#' @param ahat1 see function
#'
#' @keywords internal
#'
ahat4_func <- function(A1, A2, p, n1, n2, ahat2, ahat1){
  n <- n1 + n2 - 2
  (1 / c0_func(n)) *
    ((1 / p) * (tr(A1 + A2) ^ 4) -
       p * c1_func(n) * ahat1 -
       (p ^ 2) * c2_func(n) * (ahat1 ^ 2) * ahat2 -
       p * c3_func(n) * (ahat2 ^ 2) -
       n * (p ^ 3) * (ahat1 ^ 4))
}

#' Estimator for Frobenius norm expansion term (helper function)
#'
#' @param A1 sum of squares for group 1
#' @param A2 sum of squares for group 2
#' @param p dimension
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#' @param ahat2 see function
#' @param ahat1 see function
#'
#' @keywords internal
#'
ahat3_func <- function(A1, A2, p, n1, n2, ahat2, ahat1){
  n <- n1 + n2 - 2
  (1 / (n * ((n ^ 2) + 3 * n + 4))) *
    (((1 / p) * tr((A1 + A2) ^ 3)) -
       (3 * n * (n + 1) * p * ahat2 * ahat1) -
       (n * (p ^ 2) * (ahat1 ^ 3)))
}

#' Estimator for frobenius norm expansion term Srivastava 2014 (helper function)
#'
#' @param n sample size for groups
#' @param sample_covs covariance matrices for groups
#' @param see function
#' @param A sum of squares
#' @param overall_n overall sample size
#'
#' @keywords internal
#'
ahat2iSrivastava2014_func <- function(n, p, D, A){
  ((n - 2) * (n - 1) * tr(A ^ 2) -
     n * (n - 1) * tr(D ^ 2) +
     tr(A) ^ 2) /
    (p * n * (n - 1) * (n - 2) * (n - 3))
}

#' Estimator for frobenius norm expansion term Srivastava 2014 (helper function)
#'
#' @param ahat21 see function
#' @param ahat22 see funciton
#' @param n1 sample size group 1
#' @param n2 sample size group 2
#'
#' @keywords internal
#'
ahat2Srivastava2014_func <- function(ahat21, ahat22, n1, n2){
  ((n1 - 1) * ahat21 + (n2 - 1) * ahat22) / (n1 + n2 - 2)
}
