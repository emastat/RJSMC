#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//  @title stabilize
// Function for normalize a vector of probability using the "max" truck to avoid numerical problems
// @param weight_vec (NumericVector) non-normalized vector with log - probabilities
// @return vector with normalized probabilities
NumericVector stabilize(const NumericVector weight_vec){

  int U = weight_vec.size();

  NumericVector norm_vec(U);
  NumericVector log_norm_vec(U);


  // employe the sum-log-exp trick

  double max_value = max(weight_vec)  ;

  NumericVector scaled_vec = weight_vec - max_value ;


  double  sum_exp_scaled =0 ;

  for(int i=0 ;i < U ; i++){

    sum_exp_scaled += std::exp(scaled_vec[i]) ;

  }

  double log_sum_exp_scaled = std::log(sum_exp_scaled) ;

  //compute non-normalized weight

  for(int i=0; i<U ; i++){

    log_norm_vec[i] = weight_vec[i] - max_value - log_sum_exp_scaled ;

    norm_vec[i] = std::exp(log_norm_vec[i]) ;
  }


  return norm_vec ;

}


// @title table_cpp
// Function for computing frequency table of a given vector, for a given set of integer
// @param X (IntegerVector) Vector containing the data
// @param n (int) total number of "categories"; n=3 means the categories are ("0",1,2,3)
// @param include_zero (bool) parameters specifying if the category "0" must be taken into account
// @return frequency table
IntegerVector table_cpp(IntegerVector X, int n,bool include_zero){


  int X_length =X.size() ;

  if(X_length == 0) stop("vector of lenght 0 in table_cpp") ;

  int X_i = 0 ;

  IntegerVector count_vec(n+1);

  for(int i=0; i<X_length;i++){

    X_i = X[i];

    count_vec[ X_i ] = count_vec[ X_i ] + 1 ;

  }

  if(include_zero==false) count_vec.erase(0) ;

  return count_vec ;

}
