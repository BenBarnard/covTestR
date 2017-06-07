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
double Ishii2016Stat(List x) {
  int len = x.length();
  arma::mat pmat = x[1];
  List samplecov(len);
  double p = pmat.n_cols;
  double n = pmat.n_rows;
  double ntot = 0;
  arma::mat Apool(p, p);
  arma::mat dpool(n, n);
  arma::vec ns(len);
  List Ai(len);
  List Di(len);
  List di(len);
  List lambda(len);
  List eigendual(len);
  List lambdaTilde(len);
  List eigendualTilde(len);
  List ki(len);

 for(int i = 0; i < len; ++i){
    arma::mat mats = x[i];
    double nsi = mats.n_rows;
    ns[i] = nsi;
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
    arma::mat D = d * pow(ns - 1, -1);
    di[i] = d;
    Di[i] = D;

    arma::vec lamb;
    arma::mat eigdual;
    eig_sym(lamb, eigdual, covar);

    double lambtilde = lamb[1] - (trace(D) - lamb[1]) * pow(nsi - 2, -1);
    arma::vec htilde = pow((nsi - 1) * lambtilde, -1 / 2) * scaleddf * eigdual.col(1);
    double kii = trace(D) - lambtilde;

    ki[i] = kii;
    lambda[i] = lamb;
    eigendual[i] = eigdual;
    lambdaTilde[i] = lambtilde;
    eigendualTilde[i] = htilde;

    ntot += ns[i] - 1;
    Apool += A;
    dpool += d;
  }

  arma::mat overallcov = Apool / ntot;
  arma::mat overallD = dpool / ntot;
  arma::vec overallLambda;
  arma::mat overalleigendual;
  eig_sym(overallLambda, overalleigendual, overallcov);
  double k = trace(overallD) - overallLambda[1];

  double stat = 0;
  for(int i = 0; i < len; ++i){
    double lambda = lambdaTilde[i];
    arma::vec eigenTilde = eigendualTilde[i];
    double h = arma::as_scalar(eigenTilde.t() * overalleigendual[1]);
    double kii = ki[i];
    stat += abs(lambda * pow(overallLambda[1], -1) *
      h * kii * pow(k, -1) - 1);
  }

  return stat;
}

