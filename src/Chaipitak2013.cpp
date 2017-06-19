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
double Chaipitak2013Stat(List x) {
  int len = x.length();
  arma::vec a2i(len);
  arma::mat pmat = x[1];
  double p = pmat.n_cols;
  double ntot = 0;
  arma::mat Apool(p, p);
  Apool.fill(0);
  double ninv = 0;
  double ninv2 = 0;

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int ns = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat diag(ns, ns);
    diag.fill(0);
    diag.eye(ns, ns);
    arma::mat J(ns, ns);
    J.fill(1);
    arma::mat A = mats.t() * (diag - J / ns) * mats;
    arma::mat covar = A / (ns - 1);
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

  double f = - 4 / ntot;
  double c = - (2 * pow(ntot, 2) + 3 * ntot - 6) / (ntot * (pow(ntot, 2) + ntot + 2));
  double d =  (2 * (5 * ntot + 6)) / (ntot * (pow(ntot, 2) + ntot + 2));
  double e = - (5 * ntot + 6) / (pow(ntot, 2) * (pow(ntot, 2) + ntot + 2));

  double tau = (pow(ntot, 5) * (pow(ntot, 2) + ntot + 2)) /
    ((ntot + 1) * (ntot + 2) * (ntot + 4) * (ntot + 6) * (ntot - 1) *
    (ntot - 2) * (ntot - 3));

  double a4star = (tau / p) *
    (trace(overallcov * overallcov * overallcov * overallcov) +
    f * trace(overallcov * overallcov * overallcov) * trace(overallcov) +
    c * pow(trace(overallcov * overallcov), 2) +
    d * trace(overallcov * overallcov) * pow(trace(overallcov), 2) +
    e * pow(trace(overallcov), 4));

  double delta2 = 4 * ((2 * a4star * ninv) / (pow(a2, 2) * p) + ninv2);

  double stat = 0;
  for(int i = 0; i < len; ++i){
    stat += pow(a2i[i] / a2 - 1, 2) / delta2;
  }

  return stat;
}

