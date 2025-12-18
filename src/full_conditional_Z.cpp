#include <Rcpp.h>
#include "toolbox2.hpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' @title full_conditional_Z
//' Function for sampling the Z state for a  new non-empty segment from its full conditional and/or evaluate it at the sampled point
//' @param Z (const int&) Value of the current Z if present; if not set it 0.
//' @param K (const int&) Number of levels for the Z state
//' @param L (const double&) Length of the selected segment
//' @param keyvec (const NumericVector&) Shape parameters ruling the length of segment for the different levels of the Z state
//' @param etavec (const NumericVector&) Mean parameters ruling the length of segment for the different levels of the Z state
//' @param probvec_Z (const NumericVector&) Probability mass for the Z state distribution
//' @param sample_Z ( const bool=false) if true a new Z is sampled and the full conditional is evaluated at the new value. If false (by default) the full condititional is evaluated at the input Z
//' @return A list containting 2 elements:
//' \describe{
//' \item{Z_new}{int, the sample Z state for the selected segment}
//' \item{eval_densZ}{double, value of the  log full conditional for either the new Z sampled of the current Z inputed ( depending on sample_Z)}
//' }
//' @name full_conditional_Z
// [[Rcpp::export]]
List full_conditional_Z(const int& Z,
                        const int& K,
                        const double& L,
                        const NumericVector& keyvec,
                        const NumericVector& etavec,
                        const NumericVector& probvec_Z,
                        const bool& sample_Z=false){


  // initialize Z and other the required quantities

  int Z_new=0 ;

  double eval_densZ, log_gamma ;

  NumericVector log_densZ(K), densZ(K);

  // compute probability weights for the full conditional

  for(int i=0;i<K;i++){

    log_gamma = R::dgamma(L,keyvec[i],etavec[i]/keyvec[i],true) ;

    log_densZ[i]= std::log( probvec_Z[i] ) + log_gamma  ;

  }

  // employ the log-sum-exp trick
  densZ = stabilize(log_densZ)  ; // final probability vector


  // a new Z must be sampled - log full conditional evaluated in this point
  if(sample_Z==true){

    // sample Z

    Z_new = as<int>( Rcpp::sample( K, 1, false, densZ,true ) ) ;

    // logarithm of the evaluated the density in the sampled point
    eval_densZ = std::log(densZ[ Z_new-1 ]);

    // return a list with the new sample Z and the value of the full conditional in this value

    return( List::create(Named("Z")=Z_new,
                         Named("eval_densZ")=eval_densZ)) ;

  }else{

    // Z must NOT be sampled - log full conditional is evaluated at the current Z

    // logarithm of the evaluated the density in the current point

    eval_densZ = std::log(densZ[ Z-1 ]);



    // return a list with the current Z and the value of the full conditional in this value

    return( List::create(Named("Z")=Z,
                         Named("eval_densZ")=eval_densZ)) ;

  }

}







