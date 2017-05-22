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
double c3(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
    for(int j = i + 1; j < nrow; ++j){
      total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
        x.submat(j, 0, j, ncol - 1).t() *
        x.submat(i, 0, i, ncol - 1) *
        x.submat(j, 0, j, ncol - 1).t());
    }
  }
  return total * 2 / (nrow * (nrow - 1));
}

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
double c2(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
    for(int j = i + 1; j < nrow; ++j){
          total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
            x.submat(i, 0, i, ncol - 1).t() *
            x.submat(j, 0, j, ncol - 1) *
            x.submat(j, 0, j, ncol - 1).t());
    }
  }
  return total * 2 / (nrow * (nrow - 1));
}

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
double c1(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
        total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
          x.submat(i, 0, i, ncol - 1).t());
  }
  return total / nrow;
}
