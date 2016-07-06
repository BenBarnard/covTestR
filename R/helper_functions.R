#' Simulate Monte Carlo Samples from a Multivariate Normal Distribution
#'
#' @param meanVec Numeric vector of means for the variables.
#' @param covMat Positive-definite symmetric matrix of the variables.
#' @param maxN Numeric scalar of the maximum of the sample sizes in the simulation.
#' @param ... Other variables used in mvnorm function in the mass package.
#'
#' @return Matrix of maxN multivariate normal samples.
#'
#' @importFrom dplyr rename
#' @importFrom plyr ldply
#' @importFrom stats setNames
#' @importFrom reshape2 melt
#' @importFrom MASS mvrnorm
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' mcSamples(c(0,0,0), diag(1, 3), 10, 2)
mcSamples <- function(meanVec, covMat, maxN, pops, ...){
  replicate(pops,
            mvrnorm(n = maxN, mu = meanVec, Sigma = covMat) %>%
              melt %>%
              setNames(c('Subjects', 'Variables', 'Value')),
            simplify = FALSE) %>%
    setNames(c(1:pops)) %>%
    ldply %>%
    rename(Group = `.id`)
}

#' b hat used in method introduced by Chaipitak 2013 (helper function)
#'
#' @param ahat21 see function
#' @param ahat22 see function
#'
#' @keywords internal
#'
bhat_func <- function(ahat21, ahat22){
  ahat21 / ahat22
}

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

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
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
c1_func <- function(n){
  2 * n * (2 * (n ^ 2) + 6 * n + 9)
}

#' Sample size expansion term Srivastava 2007 (helper function)
#'
#' @param n sample size
#'
#' @keywords internal
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
c3_func <- function(n){
  n * (2 * (n ^ 2) + 5 * n + 7)
}

#' Overall Covariance Matirx (helper function)
#'
#' @param A1 sum of squares for group 1
#' @param A2 sum of squares for group 2
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#'
#' @keywords internal
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
A_func <- function(matrix){
  (t(matrix) - colMeans(matrix)) %*% t(t(matrix) - colMeans(matrix))
}

#' delta hat squared for Chaiptak 2013 (helper function)
#'
#' @param ahatstar4 see function
#' @param p dimension
#' @param ahat2 see function
#' @param n1 sample size for group 1
#' @param n2 sample size for group 2
#'
#' @keywords internal
#'
deltahat2_func <- function(ahatstar4, p , ahat2, n1, n2){
  4 * (((2 * ahatstar4) / (p * ahat2 ^ 2)) * sum(1 / (n1 - 1), 1 / (n2 - 1)) + sum(1 / ((n1 - 1) ^ 2), 1 / ((n2 - 1) ^ 2)))
}

#' Estimator for eta squared Srivastave 2007 (helper function)
#'
#' @param n sample size
#' @param p dimension
#' @param ahat4 see function
#' @param ahat2 see function
#'
#' @keywords internal
#'
etahat2i_func <- function(n, p, ahat4, ahat2){
  (4 / ((n - 1) ^ 2)) *
    (ahat2 ^ 2) *
    (1 + ((2 * (n - 1) * ahat4) / (p * ahat2 ^ 2)))
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

#' Estimator for Gamma Srivastava and Yanagihara 2010 (helper function)
#'
#' @param ahat2i see function
#' @param ahat1i see function
#'
#' @keywords internal
#'
gammahati_func <- function(ahat2i, ahat1i){
  ahat2i / (ahat1i ^ 2)
}

#' Estimator for Ksi Srivastava and Yanagihara 2010 (helper function)
#'
#' @param n sample size
#' @param p dimension
#' @param ahat1 see function
#' @param ahat2 see function
#' @param ahat3 see function
#' @param ahat4 see function
#'
#' @keywords internal
#'
ksihat2i_func <- function(n, p, ahat1, ahat2, ahat3, ahat4){
  (4 / ((n - 1) ^ 2)) * (((ahat2 ^ 2) / (ahat1 ^ 4)) +
                           ((2 * (n - 1)) / p) *
                           (((ahat2 ^ 3) / (ahat1 ^ 6)) -
                              ((2 * ahat2 * ahat3) / (ahat1 ^ 5)) +
                              (ahat4 / (ahat1 ^ 4))))
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

ahat2iSrivatava2014_func <- function(){

}

#' Turn Tidy data frame into data matrix (helper function)
#'
#' @param data tidy dataframe
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom reshape2 acast
#' @importFrom dplyr select
#'
dataDftoMatrix <- function(data){
  data %>%
    select(-Group) %>%
    acast(Subjects~Variables,
          value.var = "Value")
}

#' Trace of Matrix
#'
#' @param mat matrix
#'
#' @keywords internal
#'
tr <- function(mat){
  sum(diag(mat))
}
