#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Chaipitak2013_(List x) {

  ns <- lapply(matrix_ls, function(matrix){
    nrow(matrix)
  })

  p <- lapply(matrix_ls, function(matrix){
    ncol(matrix)
  })

  A_ls <- lapply(matrix_ls, A_func)

  sample_covs <- lapply(matrix_ls, cov)
  overall_cov <- overall_cov_func(A_ls, ns)

  ahat2i <- mapply(ahat2i_func, ns, p, sample_covs, SIMPLIFY = FALSE)
  ahat2 <- ahat2_func(ns, overall_cov, p[[1]])

  tau <- tau_func(ns)
  ahatStar4 <- ahatStar4_func(tau, p[[1]], overall_cov, ns)
  deltahat2 <- mapply(deltahat2_func, ahatStar4, p, ahat2, ns, SIMPLIFY = FALSE)
  b <- mapply(function(ahat2, ahat2i){ahat2i / ahat2}, ahat2, ahat2i, SIMPLIFY = FALSE)

  statistic <- Reduce(`+`, mapply(function(b, deltahat2){
    (b - 1) ^ 2 / deltahat2
  }, b, deltahat2, SIMPLIFY = FALSE))
  return statistic;
}

