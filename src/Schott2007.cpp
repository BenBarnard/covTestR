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
double Schott2007Stat(List x) {
  int len = x.length();
  arma::vec a2i(len);
  arma::mat pmat = x[1];
  List samplecov(len);
  double p = pmat.n_cols;
  double ntot = 0;
  arma::mat Apool(p, p);
  double ninv = 0;
  double ninv2 = 0;

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int ns = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat diag(ns, ns);
    diag.eye(ns, ns);
    arma::mat J(ns, ns);
    J.fill(1);
    arma::mat A = mats.t() * (diag - J / ns) * mats;
    arma::mat covar = A / (ns - 1);
    samplecov[i] = covar;
    p = ps;
    a2i[i] = (pow(ns - 1, 2) / (ps * (ns - 2) * (ns + 1))) *
      (trace(covar * covar) - pow(ns - 1, -1) * pow(trace(covar), 2));

    ntot += ns - 1;
    Apool += A;
    ninv += pow(ns - 1, -1);
    ninv2 += 1 / pow(ns - 1, 2);
  }

  arma::mat overallcov = Apool / ntot;
  double a2 = ((ntot * ntot) / (p * (ntot - 1) * (ntot + 2))) *
    (trace(overallcov * overallcov) - (1 / ntot) *
    trace(overallcov) * trace(overallcov));

  double theta = 0;

  for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int ns = mats.n_rows;
  theta += 2 * a2 * pow(ns - 1, -1);
  }

  double stat = 0;
  for(int i = 0; i < len; ++i){
    arma::mat sampcov = samplecov[i];
    stat += pow(a2i[i] + a2 - (2 / p) * trace(sampcov * overallcov), 2) /
      pow(theta, 2);
  }

  return stat;
}

