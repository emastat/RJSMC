#include <Rcpp.h>
#include "extract_index.hpp"
#include "toolbox2.hpp"

using namespace Rcpp;

//' @title weights_correction
//' Function for computing the weights correction based on the importance density for a segment
//' @param T_seg Vector with time stamps within the update interval
//' @param break_sample_output (const List&) List of breakpoint vectors generated with the Importance Sampling
//' @param M (const int&) number of iteration performed in the Importance Sampling
//' @param K (const int&) Number of levels of the Z state
//' @param keyvec (const NumericVector&) Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec (const NumericVector&) Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec (const NumericVector&) Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec (const NumericVector&) Mean parameters ruling the length of empty segments in different levels of F state
//' @param probvec_Z (const NumericVector&) probability mass for the Z state
//' @param probvec_F (const NumericVector&) probability mass for the F state
//' @param probvec_Q (const NumericVector&) probability mass for the Q state
//' @param alphavec (const NumericVector&) Shape parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param muvec (const NumericVector&) Mean parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param minimum_n (const int&) minimum number of observations that must be observed in a non-empty segment
//' @param W (const int&) Number of levels for the Q state
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @return correction weight (double)
//' @export
// [[Rcpp::export]]

double weigths_correction(const NumericVector& T_seg,
                           const List& break_sample_output,
                           const int& M,
                           const int& K,
                           const NumericVector& keyvec,
                           const NumericVector& etavec,
                           const NumericVector& key0vec,
                           const NumericVector& eta0vec,
                           const NumericVector& probvec_Z,
                           const NumericVector& probvec_F,
                           const NumericVector& probvec_Q,
                           const NumericVector& alphavec,
                           const NumericVector& muvec,
                           const int& W,
                           const double& P0){


  // initialize the correction Value

  double weight = 0.0 ; // final weight
  double weight_m = 0.0 ; // weight for current sample
  double weight_i = 0.0; // weight for current segment
  
  int F_levels = probvec_F.length() ;

  // initialize state value (it will be either Z or F)

  int state_new=0 ;
  int initial_state = 0 ;

  // initialize left and right bound of the segment

  double start_point = 0.0 ;
  double end_point = 0.0 ;

  // vector holding the observation counts in the segment
  NumericVector out_extract(3) ;
  double count = 0.0 ;

  // initialize the density values
  double ldens_L = 0.0 ;

  // sample if the segment is empty or not

  NumericVector prob_vec = {P0, 1-P0}; //P0-> prob the segment is empty

  //extract list with all the generated breakpoints vectors during Importance Sampling
  List Bvec_list = break_sample_output["S.samp"] ;

  // extract vector with log densities of the breakpoints vectors generated during Importance Sampling
  NumericVector IS_log_density_vec = break_sample_output["log_density"] ;

  // loop over the Importance sampling samples

  for(int m=0; m<M; m++){

    //retrieve breakpoint vector for sample m
    NumericVector Bvec_m = Bvec_list[m] ;

    //retrieve importance sampling log density for sample m
    double IS_log_density = IS_log_density_vec[m] ;


    //loop over the segments
    for(int i=0; i < (Bvec_m.length()-1); i++){

      // extract start_point and end_point
      start_point = Bvec_m[i] ;
      end_point = Bvec_m[i+1] ;

      //compute length of the segment
      double L_new = end_point - start_point ;

      // 1-> empty , 2-> non-empty
      initial_state = as<int>(Rcpp::sample(2, 1, false, prob_vec, true) ) ;

      // extract number of observations within the selected segment
      out_extract =clone( extract_index(T_seg, start_point, start_point + L_new)) ;

      // number of obs in [start_point; end_point)
      count = out_extract[0] ;

      // illegal configuration
      if (((count==0) & (initial_state==2)) | ((count>0) & (initial_state==1))){

        // reset the weight and IS_log_density to 0.0
        weight_m = 0.0 ;
        IS_log_density = 0.0;

        break;


      }else{

        //legal configuration

        if(count>0){

          // sample Z
          state_new = as<int>( Rcpp::sample( K, 1, false, probvec_Z, true ) ) ;

          // ldens_L: prior density for the length of the selected segment
          ldens_L = R::dgamma( L_new, keyvec[state_new-1],etavec[ state_new-1 ]/keyvec[state_new-1],true );

          //TIME-STAMPS POSITION: T_s1... Ts,n_s
          double ldens_tmsp =std::lgamma(count + 1) - count * std::log(L_new)  ;

          // sample Q
          int Q_new = as<int>( Rcpp::sample( W, 1, false, probvec_Q, true ) ) ;

          // ldens_count: (collapsed) log density for the count of obs inside a CLOSED segment
          double ldens_count = R::dnbinom(count,
                                          alphavec[Q_new-1],
                                          (alphavec[Q_new-1]/muvec[Q_new-1])/(alphavec[Q_new-1]/muvec[Q_new-1]+ L_new)
                                          ,true) ;

          // compute the  log weight for the current segment                      
          weight_i = ldens_tmsp + ldens_count + ldens_L ;

        }else{

          // sample F
          state_new =  as<int>( Rcpp::sample( F_levels, 1, false, probvec_F, true ) ) ;

          // ldens_L: prior density for the length of the selected segment
          ldens_L = R::dgamma( L_new, key0vec[state_new-1] ,eta0vec[state_new-1]/key0vec[state_new-1],true );

          // compute final correction value                                
          weight_i = ldens_L ;
        }
        
        // update weight_m
        weight_m = weight_m + weight_i ;

      }

    }

    // subtract the importance sampling log density from the weight for sample m
    weight_m = weight_m - IS_log_density ;

    // Cap weight_m to prevent overflow when exponentiating (similar to incremental_weight.cpp)
    if(weight_m > 700.0) {
      weight_m = 700.0;
    }

    //update final weight 
    weight = weight_m +  std::exp(weight_m) ;

  }

  // START CODE TO REMOVE
  if(NumericVector::is_na(weight)){

    stop("weights_correction.cpp --> weight is na");

  }

  if(Rcpp::traits::is_nan<REALSXP>(weight)){

    stop("weights_correction.cpp -->  weight is NaN");

  }

  if(Rcpp::traits::is_infinite<REALSXP>(weight)){

    stop("weights_correction.cpp -->  weight is infinite");

  }
  // END CODE TO REMOVE

  return weight/M ;

}
