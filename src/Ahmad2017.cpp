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
double Ahmad2017Stat(List x) {
  int len = x.length();
  arma::vec Eis(len);
  arma::vec Eijs(len * (len - 1) / 2);
  double ntot = 0;
  arma::mat pmat = x[1];
  double p = pmat.n_cols;
  arma::mat Apool(p, p);

  for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int ns = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat diag(ns, ns);
    diag.eye(ns, ns);
    arma::mat J(ns, ns);
    J.fill(1);
    arma::mat A = mats.t() * (diag - J / ns) * mats;
    ntot += ns - 1;
    Apool += A;
  }

  arma::mat pooledcov = Apool / ntot;;

  double doublesum = 0;
  double singlesum = 0;
  for(int i = 0; i < len; ++i){
    double Ei = Eis[i];
    singlesum += Ei;
    for(int j = i + 1; j < len; ++j){
      double Eij = Eijs[i];
      doublesum += Eij;
    }
  }
  double stat = a * ((len - 1) * singlesum - 2 * doublesum) / trace(pooledcov * pooledcov);

  return stat;
}
