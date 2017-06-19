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
  double ninv = 0;
  double ninv2 = 0;
  arma::vec ns(len);
  List Ai(len);
  List Di(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    ns[i] = mats.n_rows;
    int ps = mats.n_cols;
    arma::mat covar = cov(mats);
    samplecov[i] = covar;

    Ai[i] = covar * (ns[i] - 1);

    p = ps;

    arma::rowvec scaled = mean(mats);
    arma::mat scaleddf(ns[i], p);

    for(int k = 0; k < ps; ++k){
      scaleddf.col(k) = mats.col(k) - scaled(k);
    }

    arma::mat d = scaleddf * scaleddf.t();
    arma::mat D(ns[i], ns[i]);
    D.fill(0);

    for(int z = 0; z < ns[i]; ++z){
      D(z, z) = d(z, z);
    }

    Di[i] = D;
    ntot += ns[i] - 1;

    ninv += pow(ns[i] - 1, -1);
    ninv2 += pow(ns[i] - 1, -2);
  }

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

 double a2 = a2num * pow(ntot, -1);

 double theta = 2 * a2 * ninv;

arma::vec trcov(len);
  double stat = 0;
  for(int i = 0; i < len; ++i){
    arma::mat sampcovi = samplecov[i];
    for(int j = i + 1; j < len; ++j){
      arma::mat sampcovj = samplecov[j];
    stat += pow(a2i[i] + a2i[j] - (2 * pow(p, -1)) * trace(sampcovi * sampcovj), 2) *
      pow(theta, -2);
    }
  }

  return stat;
}

