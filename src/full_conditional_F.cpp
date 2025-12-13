#include <Rcpp.h>
#include "toolbox2.hpp"
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' @title full_conditional_F
//' Function for sampling the F state for a  new empty segment from its full conditional and/or evaluate it at the sampled point
//' @param F (const int&) Value of the current F if present; if not set it 0.
//' @param L (const double&) Length of the selected segment
//' @param key0vec (const NumericVector&) Shape parameters ruling the length of empty segment for the different levels of the F state
//' @param eta0vec (const NumericVector&) Mean parameters ruling the length of empty segment for the different levels of the F state
//' @param probvec_F (const NumericVector&) Probability mass for the F state distribution
//' @param sample_F ( const bool=false) if true a new F is sampled and the full conditional is evaluated at the new value. If false (by default) the full condititional is evaluated at the input F
//' @return A list containting 2 elements:
//' \describe{
//' \item{F_new}{int, the sample F state for the selected segment}
//' \item{eval_densF}{double, value of the  log full conditional for either the new F sampled of the current F inputed ( depending on sample_F)}
//' }

// [[Rcpp::export]]
List full_conditional_F(const int& F,
                        const double& L,
                        const NumericVector& key0vec,
                        const NumericVector& eta0vec,
                        const NumericVector& probvec_F,
                        const bool sample_F=false){

  int F_levels = probvec_F.size() ;

  // initialize F and other the required quantities

  int F_new=0 ;

  double eval_densF, log_gamma ;

  NumericVector log_densF(F_levels), densF(F_levels);

  // compute probability weights for the full conditional

  for(int i=0;i<F_levels;i++){

    log_gamma = R::dgamma(L,key0vec[i],eta0vec[i]/key0vec[i],true) ;

    log_densF[i]= std::log( probvec_F[i] ) + log_gamma  ;

  }

  // employ the log-sum-exp trick

  densF = stabilize(log_densF)  ; // final probability vector

  // a new F must be sampled - log full conditional evaluated in this point

  if(sample_F==true){

    // sample F

    F_new = as<int>( Rcpp::sample( F_levels, 1, false, densF,true ) ) ;

    // logarithm of the evaluated the density in the sampled point
    eval_densF = std::log(densF[ F_new-1 ]);

    // reuturn a list with the new sample F and the value of the full conditional in this value

    return( List::create(Named("F")=F_new,
                         Named("eval_densF")=eval_densF)) ;

  }else{

    // F must NOT be sampled - log full conditional is evaluated at the current F

    // logarithm of the evaluated the density in the current point

    eval_densF = std::log(densF[ F-1 ]);

    // reuturn a list with the current F and the value of the full conditional in this value

    return( List::create(Named("F")=F,
                         Named("eval_densF")=eval_densF)) ;


  }

}







