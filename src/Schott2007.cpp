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
  arma::vec ns(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    ns[i] = mats.n_rows;
    int ps = mats.n_cols;

    arma::mat covar = cov(mats);
    samplecov[i] = covar;
    p = ps;
    a2i[i] = (pow(ns[i] - 1, 2) / (ps * (ns[i] - 2) * (ns[i] + 1))) *
      (trace(covar * covar) - pow(ns[i] - 1, -1) * pow(trace(covar), 2));

    ntot += ns[i] - 1;
    ninv += pow(ns[i] - 1, -1);
    ninv2 += pow(ns[i] - 1, -2);
  }

 double a2num = 0;
 for(int i = 0; i < len; ++i){
   a2num += ns[i] * a2i[i];
 }

 double a2 = a2num * pow(ntot, -1);

  double theta = 0;

  for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int ns = mats.n_rows;
  theta += 2 * a2 * pow(ns - 1, -1);
  }

  double stat = 0;
  for(int i = 0; i < len; ++i){
    arma::mat sampcovi = samplecov[i];
    for(int j = i + 1; j < len; ++j){
      arma::mat sampcovj = samplecov[j];
    stat += pow(a2i[i] + a2i[j] - (2 / p) * trace(sampcovi * sampcovj), 2) /
      pow(theta, 2);
    }
  }

  return stat;
}

