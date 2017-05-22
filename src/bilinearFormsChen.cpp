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
double bilinearquad(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
    for(int j = i + 1; j < nrow; ++j){
      for(int k = j + 1; k < nrow; ++k){
        for(int l = k + 1; l < nrow; ++l){
    total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
      x.submat(j, 0, j, ncol - 1).t() *
      x.submat(k, 0, k, ncol - 1) *
      x.submat(l, 0, l, ncol - 1).t());
        }
      }
    }
  }
  return total * 24 / (R::gammafn(nrow + 1) / R::gammafn(nrow - 3));
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
double bilinearsquare(arma::mat x) {
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
  return total * 2 / (R::gammafn(nrow + 1) / R::gammafn(nrow - 1));
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
double bilinearcube(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
    for(int j = i + 1; j < nrow; ++j){
      for(int k = j + 1; k < nrow; ++k){
        total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
          x.submat(j, 0, j, ncol - 1).t() *
          x.submat(j, 0, j, ncol - 1) *
          x.submat(k, 0, k, ncol - 1).t());
      }
    }
  }
  return total * 6 / (R::gammafn(nrow + 1) / R::gammafn(nrow - 2));
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
double bilinearoff(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
    for(int j = i + 1; j < nrow; ++j){
        total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
          x.submat(j, 0, j, ncol - 1).t());
      }
  }
  return total * 2 / (R::gammafn(nrow + 1) / R::gammafn(nrow - 1));
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
double quadra(arma::mat x) {
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  double total = 0;
  for(int i = 0; i < nrow; ++i){
        total += arma::as_scalar(x.submat(i, 0, i, ncol - 1) *
          x.submat(i, 0, i, ncol - 1).t());
  }
  return total / nrow;
}
