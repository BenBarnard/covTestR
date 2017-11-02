#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Chaipitak2013poolStat(List x) {
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
    arma::mat covar = cov(mats);


    a2i[i] = pow(ns - 1, 2) *
      pow(ps * (ns - 2) * (ns + 1), -1)*
      (trace(covar * covar) - pow(ns - 1, -1) * pow(trace(covar), 2));

    ntot += ns - 1;
    Apool += covar * (ns - 1);
    ninv += pow(ns - 1, -1);
    ninv2 += pow(ns - 1, -2);
  }

  arma::mat pooledCov = Apool / ntot;
  double a2 = pow(ntot, 2) *
               pow(p * (ntot - 1) * (ntot + 2), -1) *
    (trace(pooledCov * pooledCov) - pow(ntot, -1) *
    pow(trace(pooledCov), 2));

  double f = - 4 * pow(ntot, -1);
  double c = - (2 * pow(ntot, 2) + 3 * ntot - 6) * pow(ntot * (pow(ntot, 2) + ntot + 2), -1);
  double d =  (2 * (5 * ntot + 6)) * pow(ntot * (pow(ntot, 2) + ntot + 2), -1);
  double e = - (5 * ntot + 6) * pow(pow(ntot, 2) * (pow(ntot, 2) + ntot + 2), -1);

  double tau = (pow(ntot, 5) * (pow(ntot, 2) + ntot + 2)) *
    pow((ntot + 1) * (ntot + 2) * (ntot + 4) * (ntot + 6) * (ntot - 1) *
    (ntot - 2) * (ntot - 3), -1);

  double a4star = tau * pow(p, -1) *
    (trace(pooledCov * pooledCov * pooledCov * pooledCov) +
    f * trace(pooledCov * pooledCov * pooledCov) * trace(pooledCov) +
    c * pow(trace(pooledCov * pooledCov), 2) +
    d * trace(pooledCov * pooledCov) * pow(trace(pooledCov), 2) +
    e * pow(trace(pooledCov), 4));

  double delta2 = 4 * ((2 * a4star * ninv) * pow(pow(a2, 2) * p, -1) + ninv2);

  double stat = 0;
  for(int i = 0; i < len; ++i){
    stat += pow(a2i[i] * pow(a2, -1) - 1, 2) * pow(delta2, -1);
  }

  return stat;
}

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
    arma::mat covar = cov(mats);

    a2i[i] = pow(ns - 1, 2) *
      pow(ps * (ns - 2) * (ns + 1), -1)*
      (trace(covar * covar) - pow(ns - 1, -1) * pow(trace(covar), 2));

    ntot += ns - 1;
    Apool += covar * (ns - 1);
    ninv += pow(ns - 1, -1);
    ninv2 += pow(ns - 1, -2);
  }

  arma::mat pooledCov = Apool / ntot;
  double a2 = pow(ntot, 2) *
    pow(p * (ntot - 1) * (ntot + 2), -1) *
    (trace(pooledCov * pooledCov) - pow(ntot, -1) *
    pow(trace(pooledCov), 2));

  double f = - 4 * pow(ntot, -1);
  double c = - (2 * pow(ntot, 2) + 3 * ntot - 6) * pow(ntot * (pow(ntot, 2) + ntot + 2), -1);
  double d =  (2 * (5 * ntot + 6)) * pow(ntot * (pow(ntot, 2) + ntot + 2), -1);
  double e = - (5 * ntot + 6) * pow(pow(ntot, 2) * (pow(ntot, 2) + ntot + 2), -1);

  double tau = (pow(ntot, 5) * (pow(ntot, 2) + ntot + 2)) *
    pow((ntot + 1) * (ntot + 2) * (ntot + 4) * (ntot + 6) * (ntot - 1) *
    (ntot - 2) * (ntot - 3), -1);

  double a4star = tau * pow(p, -1) *
    (trace(pooledCov * pooledCov * pooledCov * pooledCov) +
    f * trace(pooledCov * pooledCov * pooledCov) * trace(pooledCov) +
    c * pow(trace(pooledCov * pooledCov), 2) +
    d * trace(pooledCov * pooledCov) * pow(trace(pooledCov), 2) +
    e * pow(trace(pooledCov), 4));

  double delta2 = 4 * ((2 * a4star * ninv) * pow(pow(a2, 2) * p, -1) + ninv2);

  double stat = 0;
  for(int i = 0; i < len; ++i){
    for(int j = i + 1; j < len; ++j){
    stat += pow(a2i[i] * pow(a2i[j], -1) - 1, 2) * pow(delta2, -1);
    }
  }

  return stat;
}

