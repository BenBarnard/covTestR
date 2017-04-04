#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param n sample size
#' @param p dimension
#' @param sample_covs sample covariance matrix
#'
#' @keywords internal
#'
#'
ahat2i_func <- function(n, p, sample_covs){
  (((n - 1) ^ 2) / (p * (n - 2) * (n + 1))) *
    (tr(sample_covs %*% sample_covs) - (1 / (n - 1)) * (tr(sample_covs)) ^ 2)
}

#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param ns sample size
#' @param p dimension
#' @param overall_cov sample covariance matrix
#'
#' @keywords internal
#'
#'
ahat2_func <- function(ns, overall_cov, p){
  nall <- Reduce(`+`, lapply(ns, function(x){x - 1}))
  ((nall ^ 2) / (p * (nall - 1) * (nall + 2))) *
    (tr(overall_cov %*% overall_cov) - (1 / nall) * tr(overall_cov) ^ 2)
}

#' Estimator for expansion term in Frobenius Norm Schott 2007 (helper funciton)
#'
#' @param p dimension
#' @param sample_cov sample covariance matrix
#'
#' @keywords internal
#'
#'
ahat1_func <- ahat1i_func <- function(p, sample_cov){
  (1 / p) * tr(sample_cov)
}


#' Estimator for expansion term in Frobenius Norm Chaipitak 2013 (helper funciton)
#'
#' @param tau see function
#' @param p dimension
#' @param sample_cov sample covariance matrix
#' @param ns sample size for group 1
#'
#'
#'
#' @keywords internal
#'
#'
ahatStar4_func <- function(tau, p, sample_cov, ns){
  n <- Reduce(`+`, lapply(ns, function(x){x - 1}))
  b <- - 4 / n
  c <- - (2 * (n ^ 2) + 3 * n - 6) / (n * ((n ^ 2) + n + 2))
  d <-  (2 * (5 * n + 6)) / (n * ((n ^ 2) + n + 2))
  e <- - (5 * n + 6) / ((n ^ 2) * ((n ^ 2) + n + 2))
  (tau / p) *
    (tr(sample_cov %*% sample_cov %*% sample_cov %*% sample_cov) +
       b * tr(sample_cov %*% sample_cov %*% sample_cov) * tr(sample_cov) +
       c * (tr(sample_cov %*% sample_cov) ^ 2) +
       d * tr(sample_cov %*% sample_cov) * (tr(sample_cov) ^ 2) +
       e * tr(sample_cov) ^ 4)
}

#' Estimator for frobenius norm expansion term Srivastava 2007
#'
#' @param A sum of squares for group 1
#' @param p dimension
#' @param ns sample size for group 1
#' @param ahat2 see function
#' @param ahat1 see function
#'
#' @keywords internal
#'
#'
ahat4_func <- function(A, p, ns, ahat2, ahat1){
  nss <- Reduce(`+`, lapply(ns, function(x){x - 1}))
  As <- Reduce(`+`, A)
  (1 / c0_func(nss)) *
    ((1 / p) * (tr((As) %*% (As) %*% (As) %*% (As))) -
       p * c1_func(nss) * ahat1 -
       (p ^ 2) * c2_func(nss) * (ahat1 ^ 2) * ahat2 -
       p * c3_func(nss) * (ahat2 ^ 2) -
       nss * (p ^ 3) * (ahat1 ^ 4))
}

#' Estimator for Frobenius norm expansion term (helper function)
#'
#' @param A1 sum of squares for group 1
#' @param p dimension
#' @param ns sample size for group 1
#' @param ahat2 see function
#' @param ahat1 see function
#'
#' @keywords internal
#'
#'
ahat3_func <- function(A, p, ns, ahat2, ahat1){
  n <- Reduce(`+`, lapply(ns, function(x){x - 1}))
  As <- Reduce(`+`, A)
  (1 / (n * ((n ^ 2) + 3 * n + 4))) *
    (((1 / p) * tr((As) ^ 3)) -
       (3 * n * (n + 1) * p * ahat2 * ahat1) -
       (n * (p ^ 2) * (ahat1 ^ 3)))
}

#' Estimator for frobenius norm expansion term Srivastava 2014 (helper function)
#'
#' @param n sample size for groups
#' @param p function
#' @param A sum of squares
#' @param D overall sample size
#'
#' @keywords internal
#'
#'
ahat2iSrivastava2014_func <- function(n, p, D, A, ns){
  n_overall <- Reduce(sum, lapply(ns, function(x){nrow(x)}))
  n_overall_1 <- Reduce(sum, lapply(ns, function(x){nrow(x) - 1}))
  ((n - 2) * (n - 1) * tr(A %*% A) -
      n_overall * n_overall_1 * tr(D %*% D) +
     tr(A %*% A)) /
    (p * n * (n - 1) * (n - 2) * (n - 3))
}

#' Estimator for frobenius norm expansion term Srivastava 2014 (helper function)
#'
#' @param ahat2i see function
#' @param n sample size group 1
#'
#' @keywords internal
#'
#'
ahat2Srivastava2014_func <- function(ahat2i, n){
  Reduce(`+`, mapply(function(n, ahat2i){(n - 1) * ahat2i}, n, ahat2i, SIMPLIFY = FALSE)) /
    Reduce(`+`, lapply(n, function(x){x - 1}))
}
