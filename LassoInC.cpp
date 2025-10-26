#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  double value = std::abs(a) - lambda;
  if (value < 0) value = 0;
  return (a >= 0 ? 1 : -1) * value;
}



// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  int n = Xtilde.n_rows;
  
  arma::colvec resid = Ytilde - Xtilde * beta;
  return arma::dot(resid, resid) / (2.0 * n) + lambda * arma::accu(arma::abs(beta));
}



// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  
  if (Ytilde.n_rows != n) {
    stop("Number of rows in Xtilde and Ytilde must match");
  }
  if (lambda < 0) {
    stop("Lambda should be nonnegative.");
  }
  if (beta_start.n_rows != p) {
    stop("Length of beta_start must match number of columns in Xtilde");
  }
}  



// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
}