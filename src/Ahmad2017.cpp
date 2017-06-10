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
  arma::mat Eijs(len, len);
  double ntot = 0;
  arma::mat pmat = x[1];
  double p = pmat.n_cols;
  arma::mat Apool(p, p);
  List ni(len);

  for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    double E = 0;
    int ns = mats.n_rows;
    arma::mat diag(ns, ns);
    diag.eye(ns, ns);
    arma::mat J(ns, ns);
    J.fill(1);
    arma::mat A = mats.t() * (diag - J / ns) * mats;
    ntot += ns - 1;
    Apool += A;

    for(int k = 0; k < ns; ++k){
      for(int r = 0; r < ns; ++r){
        for(int kp = 0; kp < ns; ++kp)
          for(int rp = 0; rp < ns; ++rp){
            double e = 0;
            if (k == r | k == kp | k == rp | r == kp | r == rp | kp == rp) {
              e = 0;
            } else {
              arma::mat x = mats.t();
              e = as_scalar((x.col(k) - x.col(r)).t() * (x.col(kp) - x.col(rp)) * (x.col(k) - x.col(r)).t() * (x.col(kp) - x.col(rp)) +
                (x.col(k) - x.col(kp)).t() * (x.col(rp) - x.col(rp)) * (x.col(k) - x.col(kp)).t() * (x.col(r) - x.col(rp)) +
                (x.col(k) - x.col(rp)).t() * (x.col(kp) - x.col(r)) * (x.col(k) - x.col(rp)).t() * (x.col(kp) - x.col(r)));
            }
            E += e;
          }
      }
    }

    Eis[i] = E * pow(ns * (ns - 1) * (ns - 2) * (ns - 3), -1);
  }


  for(int i = 0; i < len; ++i){
    for(int j =  i + 1; j < len; ++j){
      arma::mat mati = x[i];
      arma::mat matj = x[j];
      arma::mat xi = mati.t();
      arma::mat xj = matj.t();
      double nis = ni[i];
      double njs = ni[j];
      double Eij = 0;
      for(int k = 0; k < nis; ++k){
        for(int r = 0; r < nis; ++r){
          for(int l = 0; l < njs; ++l)
            for(int s = 0; s < njs; ++s){
              double e = 0;
              if (k == r | l == s) {
                e = 0;
              } else {
                e = as_scalar((xi.col(k) - xi.col(r)).t() * (xj.col(l) - xj.col(s)) *
                  (xi.col(k) - xi.col(r)).t() * (xj.col(l) - xj.col(s)));
              }
              Eij += e;
            }
        }
      }

      Eijs(i, j) = Eij * pow(4 * nis * (nis - 1) * njs * (njs - 1), -1);
    }
  }

  arma::mat pooledcov = Apool / ntot;

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

  double a = 0;

  for(int i = 0; i < len; ++i){
    for(int j = 0; j + 1 < len; ++j){
      double nsi = ni[i];
      double nsj = ni[j];
      double Eij = Eijs(i, j);
        double dem = 0;
        for(int k = i; k < j + 1; ++k){
      dem += nsi * (nsi - 1) * nsj * (nsj - 1);
        }
        double asub = nsi * (nsi - 1) * nsj * (nsj - 1) * Eij * pow(dem, - 1);

      a += asub;
    }
  }


  double stat = a * ((len - 1) * singlesum - 2 * doublesum) / trace(pooledcov * pooledcov);

  return stat;
}
