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
double Schott2001Stat(List x) {
  int len = x.length();
  arma::mat pmat = x[1];
  List samplecov(len);
  double p = pmat.n_cols;
  double ntot = 0;
  arma::mat Apool(p, p);
  Apool.fill(0);
  arma::vec ns(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    int nsi = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat covar = cov(mats);
    samplecov[i] = covar;
    p = ps;
    ntot += nsi - 1;
    ns[i] = nsi;
    Apool += covar * (nsi - 1);
  }

  arma::mat pooledCov = Apool / ntot;


  double doublesum = 0;
  double singlesum = 0;
  for(int i = 0; i < len; ++i){
    double ni = ns[i] - 1;
    arma::mat Si = samplecov[i];
    singlesum += ni * pow(ntot, -1) * trace(Si * inv(pooledCov) * Si * inv(pooledCov));
    for(int j = 0; j < len; ++j){
      double nj = ns[j] - 1;
      arma::mat Sj = samplecov[j];
      doublesum += ni * nj * pow(ntot, -2) *
        trace(Si * inv(pooledCov) * Sj * inv(pooledCov));
    }
  }

  double stat = ntot * pow(2, -1) * (singlesum - doublesum);

  return stat;
}

