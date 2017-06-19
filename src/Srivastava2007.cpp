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
double Srivastava2007Stat(List x) {
  int len = x.length();
  arma::vec a2i(len);
  arma::mat pmat = x[1];
  List samplecov(len);
  double p = pmat.n_cols;
  double ntot = 0;
  arma::mat Apool(p, p);
  Apool.fill(0);
  double ninv = 0;
  double ninv2 = 0;
  arma::vec ns(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    ns[i] = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat diag(ns[i], ns[i]);
    diag.fill(0);
    diag.eye(ns[i], ns[i]);
    arma::mat J(ns[i], ns[i]);
    J.fill(1);
    arma::mat A = mats.t() * (diag - J / ns[i]) * mats;
    arma::mat covar = A / (ns[i] - 1);
    samplecov[i] = covar;
    p = ps;
    a2i[i] = (pow(ns[i] - 1, 2) / (ps * (ns[i] - 2) * (ns[i] + 1))) *
      (trace(covar * covar) - pow(ns[i] - 1, -1) * pow(trace(covar), 2));

    ntot += ns[i] - 1;
    Apool += A;
    ninv += pow(ns[i] - 1, -1);
    ninv2 += 1 / pow(ns[i] - 1, 2);
  }

  arma::mat overallcov = Apool / ntot;
  double a2 = ((ntot * ntot) / (p * (ntot - 1) * (ntot + 2))) *
    (trace(overallcov * overallcov) - (1 / ntot) *
    trace(overallcov) * trace(overallcov));

  double a1 = trace(overallcov) / p;
  double a4 = pow(ntot * (pow(ntot, 3) + 6 * pow(ntot, 2) + 21 * ntot + 18), -1) *
    (trace(Apool * Apool * Apool * Apool) / p -
                  p * (2 * ntot * (2 * pow(ntot, 2) + 6 * ntot + 9)) * a1 -
                  pow(p, 2) * (2 * ntot * (3 * ntot + 2)) * pow(a1, 2) * a2 -
                  p * (ntot * (2 * pow(ntot, 2) + 5 * ntot + 7)) * pow(a2, 2) -
                  ntot * pow(p, 3) * pow(a1, 4));

 arma::vec eta2i(len);
 double abar = 0;
 double abarnum = 0;
 double abardem = 0;


 for(int i = 0; i < len; ++i){
   eta2i[i] = 4 * pow(ns[i] - 1, -2) * pow(a2, 2) *
     (1 + 2 * (ns[i] - 1) * a4 * pow(p * pow(a2, 2), -1));
   abarnum += a2i[i] * pow(eta2i[i], -1);
   abardem += pow(eta2i[i], -1);
 }
 abar = abarnum / abardem;

  double stat = 0;
  for(int i = 0; i < len; ++i){
    stat += pow(a2i[i] - abar, 2) * pow(eta2i[i], -1);
  }

  return stat;
}

