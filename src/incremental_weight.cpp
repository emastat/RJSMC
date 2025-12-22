#include <algorithm>
#include <Rcpp.h>
#include "extract_index.hpp"
#include "open_negative_binomial.hpp"
#include "toolbox2.hpp"

using namespace Rcpp ;

//' @title Incremental_weight
//' Function for computing the log incremental weight for a particle
//' @param T_seg Vector with time stamps within the update interval
//' @param Y_seg Vector with the messages within the update interval
//' @param V (const int&) Value of the current V if present; if not set it 0.
//' @param Z (const int&) Value of the current Z if present; if not set it 0.
//' @param Q (const int&) Value of the current Q if present; if not set it 0.
//' @param F (const int&) Value of the current F if present; if not set it 0.
//' @param U (const int&) Number of levels for the V state
//' @param K (const int&) Number of levels of the Z state
//' @param W (const int&) Number of levels for the Q state
//' @param B_last Last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration
//' @param B1 (const double&) first generate changepoint, at the current SMC iteration
//' @param num_logs (const int&) Total number of unique messages that can be observed
//' @param lambdamat (const NumericMatrix&) probability weights ruling the message composition of the segments for the different levels of the V state
//' @param keyvec (const NumericVector&) Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec (const NumericVector&) Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec (const NumericVector&) Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec (const NumericVector&) Mean parameters ruling the length of empty segments in different levels of F state
//' @param alphavec (const NumericVector&) Shape parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param muvec (const NumericVector&) Mean parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param probvec_V (const NumericVector&) Probability mass for the V state distribution
//' @param probvec_Z (const NumericVector&) probability mass for the Z state
//' @param probvec_Q (const NumericVector&) Probability mass for the Q state distribution
//' @param probvec_F (const NumericVector&) probability mass for the F state
//' @param V_last_complete (const int&) V value for the last "closed" segment generated during the previous SMC iteration
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param minimum_n (const int&) minimum number of observations that must be observed in a non-empty segment
//' @param t_star, (double) estimate of the last simulated breakpoint
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the segment stops to be observed
//' @param open_segment (const bool, default is false) if false  the segment is considered totally observed. If true the segment is considered partially observed; this affects the likelihood computation
//' @export
// [[Rcpp::export]]
double incremental_weight(const NumericVector& T_seg,
                          const IntegerVector& Y_seg,
                          const int& V,
                          const int& Z,
                          const int& Q,
                          const int& F,
                          const int& U,
                          const int& K,
                          const int& W,
                          const double& B_last,
                          const double& B1,
                          const int& num_logs,
                          const NumericMatrix& lambdamat,
                          const NumericVector& keyvec,
                          const NumericVector& etavec,
                          const NumericVector& key0vec,
                          const NumericVector& eta0vec,
                          const NumericVector& alphavec,
                          const NumericVector& muvec,
                          const NumericVector& probvec_V,
                          const NumericVector& probvec_Z,
                          const NumericVector& probvec_Q,
                          const NumericVector& probvec_F,
                          const int& V_last_complete,
                          const double& P0,
                          const int& minimum_n,
                          const double& t_star,
                          const double& end_point=0.0,
                          const bool open_segment=false){


  // L_target: length of the first segment for the target density
  double L_target = B1 - B_last ;

  // L_importance: length of the first segment for the importance density
  double L_importance = B1 - t_star ;

  NumericVector out_extract_target(3);
  NumericVector out_extract_importance(3);

  if(open_segment==false){

    out_extract_target = extract_index(T_seg,B_last,B1);

    out_extract_importance = extract_index(T_seg,t_star,B1);

  }else{

    out_extract_target = extract_index(T_seg,B_last,end_point);

    out_extract_importance = extract_index(T_seg,t_star,end_point);

  }


  double count_target = out_extract_target[0];  //number of observations within (B_last, B1/end_point)
  double count_importance = out_extract_importance[0];  //number of observations within (t_star, B1/end_point)


  //inizialize the density values for the target
  double ldens_count_target = 0.0, //(collapsed) log density for the count of obs inside segment
         ldens_tmsp_target = 0.0,  // log density for the position of the time-stamps inside the segment
         ldens_event_target = 0.0, // log density  for the messages included inside the segment
         ldens_empty_target = 0.0, // log probability for observing an empty segment
         ldens_dynamics_target =0.0,// log prior density for the state V,Z,Q
         ldens_L_target =0.0 ; // log density for the length of the selected segment

  //inizialize the density values for the importance density
  double ldens_count_importance = 0.0, //(collapsed) log density for the count of obs inside segment
         ldens_tmsp_importance = 0.0,  // log density for the position of the time-stamps inside the segment
         ldens_event_importance = 0.0, // log density  for the messages included inside the segment
         ldens_empty_importance = 0.0, // log probability for observing an empty segment
         ldens_dynamics_importance =0.0,// log prior density for the state V,Z,Q
         ldens_L_importance =0.0 ; // log density for the length of the selected segment



  // 1. compute likelihood and prior for the target distribution

  if(count_target>0){ // non-empty segment

    // time-stamp positions
    ldens_tmsp_target =std::lgamma(count_target + 1) - count_target * std::log(L_target) ;

    if(open_segment==false){

      //if(count_target<minimum_n){
      // if((count_target>0) & (count_target<minimum_n)){
      //
      //   Rcout << "count_target: " << count_target << "\n" ;
      //   Rcout << "T_seg " << T_seg << "\n" ;
      //   Rcout << "B_last " << B_last << "\n" ;
      //   Rcout << "B1 " << B1 << "\n" ;
      //   Rcout << "t_star " << t_star << "\n" ;
      //   Rcout << "endpoint " << end_point << "\n" ;
      //   //count_target = minimum_n;
      //
      //   //stop("cannot be N<minimum_n in a closed segmnet:incremental_weight");
      //
      //   }

      // density for number of observation inside the closed segment
      ldens_count_target = R::dnbinom(count_target-minimum_n,
                                     alphavec[Q-1],
                                     (alphavec[Q-1]/muvec[Q-1])/(alphavec[Q-1]/muvec[Q-1]+ L_target),
                                     true) ;

    }else{

      // density for number of observations inside the open segment
     // ldens_count_target = prob_approx(count_target,
     //                                 end_point,
     //                                  alphavec[Q-1],
     //                                  muvec[Q-1],
     //                                  minimum_n,
     //                                  B_last,
     //                                 B1,
     //                                  0.0001);

      ldens_count_target = R::dnbinom(0,
                                      alphavec[Q-1],
                                              (alphavec[Q-1]/muvec[Q-1])/(alphavec[Q-1]/muvec[Q-1]+ L_target),
                                              true) ;

    }

    //messages inside the segment
    IntegerVector log_vector_target;
    if(out_extract_target[2] < 0){
      log_vector_target = IntegerVector(0);
    } else {
      int end_idx = (int)out_extract_target[2];
      int start_idx = (int)out_extract_target[1];
      if(end_idx >= Y_seg.size()) end_idx = Y_seg.size() - 1;
      if(start_idx < 0) start_idx = 0;
      if(start_idx > end_idx){
        log_vector_target = IntegerVector(0);
      } else {
        log_vector_target = Y_seg[Range(start_idx, end_idx)];
      }
    }

    // frequency table for the messages
    IntegerVector dj_vec_target = table_cpp(log_vector_target,num_logs,false) ;

    // density for the messages inside the segment
    for(int j=0 ; j < num_logs ; j++){

      ldens_event_target = ldens_event_target + std::log(lambdamat(V-1,j) ) * (double)dj_vec_target[j] ;

    }


    // density for the state V-Q-Z
    if(V_last_complete==0){

      ldens_dynamics_target = std::log(probvec_V[V-1]) + std::log(probvec_Z[Z-1]) + std::log(probvec_Q[Q-1]) ;

    }else{

      ldens_dynamics_target = std::log(1-P0) + std::log(probvec_V[V-1]) + std::log(probvec_Z[Z-1]) + std::log(probvec_Q[Q-1]) ;

    }

    // length of the segment
    ldens_L_target = R::dgamma( L_target, keyvec[Z-1],etavec[ Z-1 ]/keyvec[Z-1],true );

  }else{ // empty segment

    ldens_empty_target = std::log(probvec_F[F-1]) ;

    ldens_dynamics_target = std::log(P0) ;

    ldens_L_target = R::dgamma( L_target, key0vec[F-1] ,eta0vec[ F-1 ]/key0vec[F-1],true );

  }

  // 2. compute likelihood and prior for the importance distribution

  if(count_importance>0){ // non-empty segment

    // time-stamp positions
    ldens_tmsp_importance =std::lgamma(count_importance + 1) - count_importance * std::log(L_importance) ;

    if(open_segment==false){

      if(count_importance<minimum_n){stop("cannote be N<minimum_n in a closed segmnet:196");}

      // number of observation inside the segment
      ldens_count_importance = R::dnbinom(count_importance-minimum_n,
                                          alphavec[Q-1],
                                          (alphavec[Q-1]/muvec[Q-1])/(alphavec[Q-1]/muvec[Q-1]+ L_importance),
                                              true) ;

    }else{

      // number of observations inside the segment
      // ldens_count_importance = prob_approx(count_importance,
      //                                      end_point,
      //                                      alphavec[Q-1],
      //                                      muvec[Q-1],
      //                                      minimum_n,
      //                                      t_star,
      //                                      B1,
      //                                      0.0001);

      ldens_count_importance = R::dnbinom(count_importance-minimum_n,
                                          alphavec[Q-1],
                                                  (alphavec[Q-1]/muvec[Q-1])/(alphavec[Q-1]/muvec[Q-1]+ L_importance),
                                                  true) ;

    }

    //messages inside the segment
    IntegerVector log_vector_importance;
    if(out_extract_importance[2] < 0){
      log_vector_importance = IntegerVector(0);
    } else {
      int end_idx = (int)out_extract_importance[2];
      int start_idx = (int)out_extract_importance[1];
      if(end_idx >= Y_seg.size()) end_idx = Y_seg.size() - 1;
      if(start_idx < 0) start_idx = 0;
      if(start_idx > end_idx){
        log_vector_importance = IntegerVector(0);
      } else {
        log_vector_importance = Y_seg[Range(start_idx, end_idx)];
      }
    }

    // frequency table for the messages
    IntegerVector dj_vec_importance = table_cpp(log_vector_importance,num_logs,false) ;

    for(int j=0 ; j < num_logs ; j++){

      ldens_event_importance  = ldens_event_importance + std::log(lambdamat(V-1,j) ) * (double)dj_vec_importance[j] ;

    }

    // state V, Z, Q, F: note that we are assuming that V_last_complete is always >0 for the importance density
    ldens_dynamics_importance = std::log(1-P0) + std::log(probvec_V[V-1]) + std::log(probvec_Z[Z-1]) + std::log(probvec_Q[Q-1]) ;


    // length of the segment
    ldens_L_importance = R::dgamma( L_importance, keyvec[Z-1],etavec[ Z-1 ]/keyvec[Z-1],true );

  }else{ // empty segment

    ldens_empty_importance = std::log(probvec_F[F-1]) ;

    ldens_dynamics_importance = std::log(P0) ;

    ldens_L_importance = R::dgamma( L_importance, key0vec[F-1] ,eta0vec[ F-1 ]/key0vec[F-1],true );

  }

  // initialize output variable
  double out = 0.000001;

  if((count_target>0) & (count_target<minimum_n) & (open_segment==false)){

    return out;

  }else{

    out = ldens_tmsp_target - ldens_tmsp_importance
               + ldens_count_target - ldens_count_importance
               + ldens_event_target - ldens_event_importance
               + ldens_dynamics_target - ldens_dynamics_importance
               + ldens_empty_target - ldens_empty_importance
               + ldens_L_target - ldens_L_importance ;

    if(out>700){out = 700;}


    // START CODE TO DELETE
    if(R_isnancpp(out) || Rcpp::traits::is_infinite<REALSXP>(std::exp(out))){
      Rcout << "ldens_tmsp_target: " << ldens_tmsp_target << "\n" ;
      Rcout << "ldens_tmsp_importance: " << ldens_tmsp_importance << "\n" ;
      Rcout << "ldens_count_target: " << ldens_count_target << "\n" ;
      Rcout << "ldens_count_importance: " << ldens_count_importance << "\n" ;
      Rcout << "ldens_event_target: " << ldens_event_target << "\n" ;
      Rcout << "ldens_event_importance: " << ldens_event_importance << "\n" ;
      Rcout << "ldens_dynamics_target: " << ldens_dynamics_target << "\n" ;
      Rcout << "ldens_dynamics_importance: " << ldens_dynamics_importance << "\n" ;
      Rcout << "ldens_empty_target: " << ldens_empty_target << "\n" ;
      Rcout << "ldens_empty_importance: " << ldens_empty_importance << "\n" ;
      Rcout << "ldens_L_target: " << ldens_L_target << "\n" ;
      Rcout << "ldens_L_importance: " << ldens_L_importance << "\n" ;
      Rcout << "count_target: " << count_target << "\n" ;
      Rcout << "count_importance: " << count_importance << "\n" ;
      Rcout << "F:  " << F << "\n" ;
      Rcout << "B1: : " << B1 << "\n" ;
      Rcout << "B_last: : " << B_last << "\n" ;
      Rcout << "t_star: " << t_star << "\n" ;
      Rcout << "out: " << out << "\n" ;

      stop("incremental_weigth.cpp--> the weight is Nan: ");

      }
    // START CODE TO DELETE

    return out ;
  }

  // end of the function
}








