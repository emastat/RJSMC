#include <Rcpp.h>
#include "toolbox2.hpp"
#include "extract_index.hpp"
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
//' @title full_conditional_V
//' Function for sampling the V  state for a  new non-empty segment from its full conditional and/or evaluate it at the sampled point
//' @param V (const int&) Value of the current V if present; if not set it 0.
//' @param U (const int&) Number of levels for the V state
//' @param num_logs (const int&) Total number of unique messages that can be observed
//' @param T_seg (const NumericVector&) Vector with time stamps within the update interval
//' @param Y_seg (const IntegerVector&) Vector with the messages within the update interval
//' @param LB (const double&) lower bound of the selected segment
//' @param UB (const double&) upper bound of the selected segment
//' @param lambdamat (const NumericMatrix&) probability weights ruling the message composition of the segments for the different levels of the V state
//' @param probvec_V (const NumericVector&) Probability mass for the V state distribution
//' @param V_left  V value for the previous segment
//' @param open_segment (const bool, default is false) if false  the segment is considered totally observed. If true the semgent is considered partially observed; this affects the likelihood computation
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the stops to be observed
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param sample_V ( const bool=false) if true a new V is sampled and the full conditional is evaluated at the new value. If false (by default) the full condititional is evaluated at the input V
//' @return A list containting 2 elements:
//' \describe{
//' \item{V_new}{int, the sample Vstate for the selected segment or the one inputed}
//' \item{eval_densV}{double, value of the log full conditional for either the new sampled V or the current Q inputed ( depending on sample_V)}
//' }

// [[Rcpp::export]]
List full_conditional_V(const int& V,
                        const int& U,
                        const int& num_logs,
                        const NumericVector& T_seg,
                        const IntegerVector& Y_seg,
                        const double& LB,
                        const double& UB,
                        const NumericMatrix& lambdamat,
                        const NumericVector& probvec_V,
                        int V_left,
                        const double& P0,
                        const double& end_point=0.0,
                        const bool open_segment=false,
                        const bool sample_V=false){


  // initialize V and other the required quantities

  int V_new=0 ;

  double eval_densV, weight, prior_V ;

  NumericVector log_densV(U), densV(U);

  NumericVector prod_dj_logLL(num_logs) ;


  //initialize out_extract: vector with indexes of the observations contained in the segment
  NumericVector out_extract(3) ;

  // compute log_vector: vector with messages within the selected segment

  if(open_segment==false){

    // observations found in the whole segment
    out_extract =clone( extract_index(T_seg,LB,UB)) ;

  }else{

    // observations found in the partial segment
    out_extract =clone( extract_index(T_seg,LB,end_point)) ;
  }

  IntegerVector log_vector = Y_seg[Range(out_extract[1],out_extract[2])] ;



  // frequency table for the messages
  IntegerVector dj_vec = table_cpp(log_vector,num_logs,false) ;


  // compute probability weights for the full conditional

  for( int i=0; i<U ; i++){

    for(int j=0 ; j<num_logs ; j++){

      prod_dj_logLL[j] =  std::log( lambdamat(i,j) ) * (double)dj_vec[j] ;

    }

    weight = sum(prod_dj_logLL) ;

    if( V_left==0){

      prior_V = probvec_V[i] ;

    }else{

      prior_V =(1-P0) * probvec_V[i] ;

    }

    log_densV[i] = std::log(prior_V) + weight  ;

  }


  // employ the log-sum-exp trick

  densV = stabilize(log_densV) ; // final probability vector


  if(sample_V==true){

    // a new Q must be sampled - log full conditional evaluated in this point

    // SAMPLE V
    V_new = as<int>( Rcpp::sample( U, 1, false, densV, true ) ) ;

    // logarithm of the evaluated the density in the sampled point
    eval_densV = std::log(densV[V_new - 1 ]) ;

    // return a list with the new sample Q and the value of the full conditional in this value

    return( List::create(Named("V")= V_new ,
                         Named("eval_densV")=eval_densV)) ;


  }else{

    // V must NOT be sampled - log full conditional is evaluated at the current V

    // logarithm of the evaluated the density in the current point

    eval_densV = std::log(densV[ V-1 ]);

    // return a list with the current Q and the value of the full conditional in this value

    return( List::create(Named("V")=V,
                         Named("eval_densV")=eval_densV)) ;

  }


}
