// Most codes of C++ functions are based on Dirk Eddelbuettel's 
// 'RcppArmadillo' examples. RISC package rebuilds them to support 
// the calculation of sparse matrices.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::sp_mat sqrt_sp(arma::sp_mat X) {
  return sqrt(X);
}

// [[Rcpp::export]]
arma::mat cent_sp_d(arma::sp_mat X) {
  arma::mat Y = conv_to<arma::mat>::from(X);
  int c = Y.n_cols;
  rowvec meanx(c);
  meanx = mean(Y, 0);
  for(int j=0; j<c; j++){
    Y.col(j) = Y.col(j) - meanx(j);
  }
  return Y;
}

// [[Rcpp::export]]
arma::sp_mat multiply_sp_sp(arma::sp_mat X, arma::sp_mat Y) {
  return X * Y;
}

// [[Rcpp::export]]
arma::mat multiply_sp_d(arma::sp_mat X, arma::sp_mat Y) {
  arma::mat result(X * Y);
  return result;
}

// [[Rcpp::export]]
arma::mat multiply_d_d(arma::mat X, arma::mat Y) {
  return X * Y;
}

// [[Rcpp::export]]
arma::sp_mat multiply_sp_d_sp(arma::sp_mat X, arma::mat Y) {
  arma::sp_mat result(X * Y);
  return result;
}

// [[Rcpp::export]]
arma::vec multiply_sp_d_v(arma::sp_mat X, arma::mat Y) {
  arma::mat result(X * Y);
  arma::vec Z = result.as_col();
  return Z;
}

// [[Rcpp::export]]
arma::sp_mat crossprod_sp_sp(arma::sp_mat X, arma::sp_mat Y) {
  return trans(X) * Y;
}

// [[Rcpp::export]]
arma::mat crossprod_sp_d(arma::sp_mat X, arma::sp_mat Y) {
  arma::mat result(trans(X) * Y);
  return result;
}

// [[Rcpp::export]]
arma::mat crossprod_d_d(arma::mat X, arma::mat Y) {
  return trans(X) * Y;
}

// [[Rcpp::export]]
arma::sp_mat tcrossprod_sp_sp(arma::sp_mat X, arma::sp_mat Y) {
  return X * trans(Y);
}

// [[Rcpp::export]]
arma::mat tcrossprod_d_d(arma::mat X, arma::mat Y) {
  return X * trans(Y);
}

// [[Rcpp::export]]
arma::rowvec winsorize_(arma::rowvec x, double y) {
  arma::uvec id = find(x >= y);
  x.elem(id).fill(y);
  return x;
}

// [[Rcpp::export]]
arma::colvec lm_coef(arma::mat X, arma::colvec y) {
  arma::colvec coef = arma::solve(X, y);
  return coef;
}

// [[Rcpp::export]]
Rcpp::List lm_(arma::mat X, arma::colvec y) {
  
  int n = X.n_rows, k = X.n_cols;
  arma::colvec coef = arma::solve(X, y);
  arma::colvec res = y - X * coef;
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0) / (n - k);
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X) * X)));  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr") = std_err,
                            Rcpp::Named("df.residual") = n - k);
}


