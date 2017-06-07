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
double Srivastava2014Stat(List x) {
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
  List Ai(len);
  List Di(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    ns[i] = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat diag(ns[i], ns[i]);
    diag.eye(ns[i], ns[i]);
    arma::mat J(ns[i], ns[i]);
    J.fill(1);
    arma::vec j(ns[i]);
    j.fill(1);
    arma::mat A = mats.t() * (diag - J / ns[i]) * mats;
    Ai[i] = A;
    arma::mat covar = A / (ns[i] - 1);
    samplecov[i] = covar;
    p = ps;
    arma::mat scaled = mats.t() * j / 50;
    arma::mat scaleddf(ns[i], ns[i]);
    for(int i = 0; i < ps; ++i){
      scaleddf.col(i) = mats.col(i) - scaled(i);
    }
    arma::mat d = scaleddf * scaleddf.t();
    arma::mat D(ps, ps);
    for(int i = 0; i < ps; ++i){
      D(i, i) = d(i, i);
    }
    Di[i] = D;
    ntot += ns[i] - 1;
    Apool += A;
    ninv += pow(ns[i] - 1, -1);
    ninv2 += 1 / pow(ns[i] - 1, 2);
  }

  arma::mat overallcov = Apool / ntot;

 for(int i = 0; i < len; ++i){
   arma::mat Ais = Ai[i];
   arma::mat D = Di[i];
   a2i[i] = pow(p * ns[i] * (ns[i] - 1) * (ns[i] - 2) * (ns[i] - 3), -1) *
     ((ns[i] - 2) * (ns[i] - 1) * trace(Ais * Ais) -
     (ntot + len) * ntot * trace(D * D) +
     trace(Ais * Ais));
 }


 double a2num = 0;
 for(int i = 0; i < len; ++i){
   a2num += ns[i] * a2i[i];
 }

 double a2 = a2num / ntot;

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

