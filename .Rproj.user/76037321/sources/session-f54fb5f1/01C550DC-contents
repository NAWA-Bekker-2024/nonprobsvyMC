#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double cfgauss_validation_cpp(NumericVector beta, List data) {
  NumericVector y = data["y"];
  NumericMatrix x = data["x"];
  NumericMatrix setm = data["setm"];
  NumericMatrix w = data["w"];
  int n = y.size();
  int px = x.ncol();
  int ms = setm.nrow();
  int pm = setm.ncol();

  double sigma2 = -2 * pow(beta[px+pm+1], 2);
  double lnsqrt2pi = log(sqrt(2*M_PI) * beta[px+pm+1]);

  NumericVector eta = beta[0] + x * beta[Range(1, px)];
  NumericMatrix eta_mat(n, ms);
  for(int j = 0; j < ms; j++) {
    eta_mat(_,j) = eta + setm(j,_) * beta[Range(px+1, px+pm)];
  }

  NumericVector tmp(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ms; j++) {
      tmp[i] += exp(pow(y[i] - eta_mat(i,j), 2) / sigma2) * w(i,j);
    }
  }

  double ret = sum(log(tmp)) - n*lnsqrt2pi;
  return -ret;
}

// [[Rcpp::export]]
NumericVector cggauss_validation_cpp(NumericVector beta, List data) {
  NumericVector y = data["y"];
  NumericMatrix x = data["x"];
  NumericMatrix setm = data["setm"];
  NumericMatrix w = data["w"];
  int n = y.size();
  int px = x.ncol();
  int ms = setm.nrow();
  int pm = setm.ncol();

  double sigma2a = pow(beta[px+pm+1], 2);
  double sigma2 = -2 * sigma2a;

  NumericVector eta = beta[0] + x * beta[Range(1, px)];
  NumericMatrix eta_mat(n, ms);
  for(int j = 0; j < ms; j++) {
    eta_mat(_,j) = eta + setm(j,_) * beta[Range(px+1, px+pm)];
  }

  NumericMatrix predicts = rep(y, ms) - eta_mat;
  NumericMatrix exppredicts = exp(pow(predicts, 2) / sigma2) * w;

  NumericVector tempsumm = rowSums(exppredicts);
  NumericVector tmp = rowSums(predicts * exppredicts) / tempsumm;

  NumericVector ret(px+pm+2);
  ret[0] = sum(tmp);
  for(int k = 1; k <= px; k++) {
    ret[k] = sum(tmp * x(_,k-1));
  }
  for(int k = px+1; k <= px+pm; k++) {
    ret[k] = sum(rowSums(predicts * exppredicts * rep(setm(_,k-px-1), each=n)) / tempsumm);
  }
  ret[px+pm+1] = sum((rowSums(pow(predicts, 2) * exppredicts) / sigma2a - 1) / tempsumm);

  for(int k = 0; k <= px+pm; k++) {
    ret[k] /= sigma2a;
  }
  ret[px+pm+1] /= beta[px+pm+1];

  return -ret;
}

// [[Rcpp::export]]
double cflogit_validation_cpp(NumericVector beta, List data) {
  NumericVector y = data["y"];
  NumericMatrix x = data["x"];
  NumericMatrix setm = data["setm"];
  NumericMatrix w = data["w"];
  int n = y.size();
  int px = x.ncol();
  int ms = setm.nrow();
  int pm = setm.ncol();

  NumericVector eta = beta[0] + x * beta[Range(1, px)];
  NumericMatrix eta_mat(n, ms);
  for(int j = 0; j < ms; j++) {
    eta_mat(_,j) = eta + setm(j,_) * beta[Range(px+1, px+pm)];
  }

  NumericMatrix p = 1 / (1 + exp(-eta_mat));

  NumericVector tmp(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ms; j++) {
      tmp[i] += (y[i] == 1 ? p(i,j) : 1-p(i,j)) * w(i,j);
    }
  }

  double ret = sum(log(tmp));
  return -ret;
}

// [[Rcpp::export]]
NumericVector cglogit_validation_cpp(NumericVector beta, List data) {
  NumericVector y = data["y"];
  NumericMatrix x = data["x"];
  NumericMatrix setm = data["setm"];
  NumericMatrix w = data["w"];
  int n = y.size();
  int px = x.ncol();
  int ms = setm.nrow();
  int pm = setm.ncol();

  NumericVector eta = beta[0] + x * beta[Range(1, px)];
  NumericMatrix eta_mat(n, ms);
  for(int j = 0; j < ms; j++) {
    eta_mat(_,j) = eta + setm(j,_) * beta[Range(px+1, px+pm)];
  }

  NumericMatrix p = 1 / (1 + exp(-eta_mat));

  NumericVector lik(n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < ms; j++) {
      lik[i] += (y[i] == 1 ? p(i,j) : 1-p(i,j)) * w(i,j);
    }
  }

  NumericMatrix tmp = (y == 1 ? 1.0 : -1.0) * p * (1-p) * w;

  NumericVector ret(px+pm+1);
  ret[0] = sum(rowSums(tmp) / lik);
  for(int k = 1; k <= px; k++) {
    ret[k] = sum((rowSums(tmp) / lik) * x(_,k-1));
  }
  for(int k = px+1; k <= px+pm; k++) {
    ret[k] = sum(rowSums(tmp / rep(lik, ms) * rep(setm(_,k-px-1), each=n)));
  }

  return -ret;
}
