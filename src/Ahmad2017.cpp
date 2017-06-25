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
  double ntot = 0;
  arma::mat pmat = x[0];
  double p = pmat.n_cols;
  arma::mat Apool(p, p);
  Apool.fill(0);
  double Ei = 0;
  double Eij = 0;
  double ninv = 0;
  double nijinv = 0;

  for(int i = 0; i < len; ++i){
    arma::mat mati = x[i];
    double ni = mati.n_rows;
    arma::mat covar = cov(mati);
    ntot += ni - 1;
    Apool += covar * (ni - 1);
    double E = 0;

    for(int k = 0; k < ni; ++k){
      for(int r = k + 1; r < ni; ++r){
        for(int kp = r + 1; kp < ni; ++kp){
          for(int rp = kp + 1; rp < ni; ++rp){
            E += as_scalar(pow((mati.row(k).t() - mati.row(r).t()).t() * (mati.row(kp).t() - mati.row(rp).t()), 2) +
            pow((mati.row(k).t() - mati.row(kp).t()).t() * (mati.row(r).t() - mati.row(rp).t()), 2) +
            pow((mati.row(k).t() - mati.row(rp).t()).t() * (mati.row(kp).t() - mati.row(r).t()), 2));;
            }
          }
        }
      }

    Ei += 2 * E * pow(ni * (ni - 1) * (ni - 2) * (ni - 3), -1);

    ninv += pow(ni, -2);

    for(int j =  i + 1; j < len; ++j){
      arma::mat matj = x[j];
      double nj = matj.n_rows;
      double Eijs = 0;
      for(int k = 0; k < ni; ++k){
        for(int r = k + 1; r < ni; ++r){
          for(int l = 0; l < nj; ++l) {
            for(int s = l + 1; s < nj; ++s){
              Eijs += as_scalar(pow((mati.row(k).t() - mati.row(r).t()).t() *
              (mati.row(l).t() - mati.row(s).t()), 2));;
            }
          }
        }
      }

      Eij += Eijs * pow(ni * (ni - 1) * nj * (nj - 1), -1);

      nijinv += 2 * pow(ni * nj, -1);
    }
  }

  arma::mat pooledCov = Apool / ntot;

  double stat = (len - 1) * Ei - 2 * Eij * pow(4 * (pow(len - 1, 2) * ninv + nijinv), -0.5) * pow(trace(pooledCov * pooledCov), -1);

  return stat;
}
