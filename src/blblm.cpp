#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' fastlm
//' @param X mat
//' @param y vector
//' @param wi vector
//' @export
// [[Rcpp::export]]
List fastlm(const arma::mat& X, const arma::colvec& y, const arma::vec& wi) {
  int n = X.n_rows, k = X.n_cols;

  arma::mat weight = diagmat(wi); //get the weight vector

  arma::colvec coef = arma::solve(((arma::trans(X))*weight*X) , ((arma::trans(X)*weight)*y));
  // fit model y ~ X

  arma::colvec res  = y - X*coef; // residuals

  double weight_sum = accu(wi); //calculate weight sum

  double num = accu( weight * (arma::pow(res,2.0)) ); //numerator of sigma
  double denom = weight_sum - k; //denominator of sigma
  double s =  sqrt(num / denom); //weighted sigma

  //outputs of the variables
  return List::create(Named("coefficients") = coef,
                      Named("weight sum")   = weight_sum,
                      Named("sigma")         = s,
                      Named("df.residual")  = n - k,
                      Named("rank")         = k,
                      Named("response")     = y,
                      Named("Res")  = res,
                      Named("weight") = wi);
}