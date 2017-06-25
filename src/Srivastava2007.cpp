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
  arma::vec ns(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int n = mats.n_rows;
    ns[i] = n;
    int ps = mats.n_cols;
    arma::mat covar = cov(mats);
    samplecov[i] = covar;


    a2i[i] = pow(n - 1, 2) *
      pow(ps * (n - 2) * (n + 1), -1)*
      (trace(covar * covar) - pow(n - 1, -1) * pow(trace(covar), 2));

    ntot += n - 1;
    Apool += covar * (n - 1);
  }

 arma::mat pooledCov = Apool / ntot;
 double a2 = pow(ntot, 2) *
   pow(p * (ntot - 1) * (ntot + 2), -1) *
   (trace(pooledCov * pooledCov) - pow(ntot, -1) *
   pow(trace(pooledCov), 2));

  double a1 = trace(pooledCov) * pow(p, -1);

  double c0 = pow(ntot, 4) + 6 * pow(ntot, 3) + 21 * pow(ntot, 2) + 18 * ntot;
  double c1 = 4 * pow(ntot, 3) + 12 * pow(ntot, 2) + 18 * ntot;
  double c2 = 6 * pow(ntot, 2) + 4 * ntot;
  double c3 = 2 * pow(ntot, 3) + 5 * pow(ntot, 2) + 7 * ntot;

  double a4 = pow(c0, -1) *
    (pow(p, -1) * trace(Apool * Apool * Apool * Apool) -
    p * c1 * a1 -
    pow(p, 2) * c2 * pow(a1, 2) * a2 -
    p * c3 * pow(a2, 2) -
    ntot * pow(p, 3) * pow(a1, 4));

 arma::vec eta2i(len);
 double abarnum = 0;
 double abardem = 0;

 for(int i = 0; i < len; ++i){
   eta2i[i] = 4 * pow(ns[i] - 1, -2) * pow(a2, 2) *
     (1 + 2 * (ns[i] - 1) * a4 *
     pow(p * pow(a2, 2), -1));

   abarnum += a2i[i] * pow(eta2i[i], -1);
   abardem += pow(eta2i[i], -1);
 }

 double abar = abarnum / abardem;

  double stat = 0;
  for(int i = 0; i < len; ++i){
    stat += pow(a2i[i] - abar, 2) * pow(eta2i[i], -1);
  }

  return stat;
}

