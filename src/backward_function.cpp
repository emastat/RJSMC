#include <algorithm>
#include <Rcpp.h>
#include "segment_box.hpp"
#include "extract_index.hpp"
#include "updating_vec.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @title backward_function
//' This function takes in the current state of the MCMC and employ a forward step where a split of an interval is proposed
//' all the input parameters described below refers to their  ost recent values.
//' @param S (int). The total number of interarrivals (intervals).
//' @param Tvec (const NumericVector). Vector of the observed time stamps in the update interval
//' @param Yvec (const IntegerVector). Vector of the observed messages in the update interval
//' @param U (const int&) Number of levels for the V state
//' @param K (const int&) Number of levels of the Z state
//' @param W (const int&) Number of levels for the Q state
//' @param probvec_V (const NumericVector&) Probability mass for the V state distribution
//' @param probvec_Z (const NumericVector&) probability mass for the Z state
//' @param probvec_Q (const NumericVector&) Probability mass for the Q state distribution
//' @param probvec_F (const NumericVector&) probability mass for the F state
//' @param lambdamat (const NumericMatrix&) probability weights ruling the message composition of the segments for the different levels of the V state
//' @param keyvec (const NumericVector&) Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec (const NumericVector&) Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec (const NumericVector&) Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec (const NumericVector&) Mean parameters ruling the length of empty segments in different levels of F state
//' @param alphavec (const NumericVector&) Shape parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param muvec (const NumericVector&) Mean parameters ruling the rate of occurrence of segment for the different levels of the Q state.
//' @param num_logs (const int). Total numer of different logs (levels of the factor )
//' @param Bvec (NumericVector). Current vector of changepoints
//' @param Zvec (IntegerVector). Current vector with Z state of each  segment
//' @param Qvec (IntegerVector). Current vector with Q state of each  segment
//' @param Vvec (IntegerVector). Current vector with V state of each  segment
//' @param Fvec (IntegerVector). Current vector with F state of each  segment
//' @param Jss1 (const double). Probability to propose a forward move
//' @param Js1s (const double). Probability to propose a backward move
//' @param P0 (const double&) probability of observing a non-empty segment after a empty one
//' @param minimum_n (const int&) minimum number of observations in a non-empty segment
//' @param start_point (double, default is 0.0) start of the update interval
//' @param end_point (double, default is 0.0) end of the update interval; to be used for partially observed segments


void       backward_function(      int& S,
                              const NumericVector& Tvec,
                              const IntegerVector& Yvec,
                              const int& U,
                              const int& K,
                              const int& W,
                              const NumericVector& probvec_V,
                              const NumericVector& probvec_Z,
                              const NumericVector& probvec_Q,
                              const NumericVector& probvec_F,
                              const NumericMatrix& lambdamat,
                              const NumericVector& keyvec,
                              const NumericVector& etavec,
                              const NumericVector& key0vec,
                              const NumericVector& eta0vec,
                              const NumericVector& alphavec,
                              const NumericVector& muvec,
                              const int& num_logs,
                                    NumericVector& Bvec,
                                    IntegerVector& Zvec,
                                    IntegerVector& Qvec,
                                    IntegerVector& Vvec,
                                    IntegerVector& Fvec,
                              const double& Jss1,
                              const double& Js1s,
                              const double& P0,
                              const int&  minimum_n,
                              const double& start_point=0.0,
                              const double& end_point=0.0){

  //Rcout << "start backward " << "\n" ;


  //initialize quantities needed for performing the backward move

  int BR ; // index (from 1 to S-1) of the breakpoint to remove

  double LB=0,  // lower bound of the selected "double" segment
         UB=0,  // upper bound of the selected "double" segment
         Td=0,  // the breakpoint that will be possibly removed
         L_12p=0; // length of the new "merged" segment

  int Z_1 =0, // Z state of the current "left" segment
      V_1 =0, // V state of the current "left" segment
      Q_1 =0, // Q state of the current "left" segment
      F_1 =0, // F state of the current "left" segment
      Z_2 =0, // Z state of the current "right" segment
      V_2 =0, // V state of the current "right" segment
      Q_2 =0, // Q state of the current "right" segment
      F_2 =0, // F state of the current "right" segment
      V_left=0; //V  state of the segment before the current one

  // open_segment: is TRUE if the merge includes the last segment
  bool open_segment = false ;


  if(S>2){

    BR = as<int>( Rcpp::sample(S-1, 1)) ;

  }else{ BR=1 ;}


  LB = Bvec[BR-1]  ;

  UB = Bvec[BR+1] ;

  Td = Bvec[BR] ;


  V_1 = Vvec[BR-1]  ;
  V_2 = Vvec[BR] ;

  Z_1 = Zvec[BR-1]  ;
  Z_2 = Zvec[BR]   ;

  Q_1 = Qvec[BR-1]  ;
  Q_2 = Qvec[BR]   ;

  F_1 = Fvec[BR-1]  ;
  F_2 = Fvec[BR]   ;


  if(BR>1){

    V_left = Vvec[ BR - 2 ] ;

  }else{

    V_left = -1 ;
  }


  if(BR==S-1){

    open_segment =true ;
  }


  L_12p = UB - LB ;

  // 1. sample states for the  new "merged" segment: V_12p, Z_12p, Q_12p, F_12p
  // 2. evaluate the proposal distribution at the new sampled values
  // 3. evaluate the posterior of the new "merged" segment

  List out_segment_12p = segment_box(Tvec,
                                     Yvec,
                                     0,
                                     0,
                                     0,
                                     0,
                                     U,
                                     K,
                                     W,
                                     LB,
                                     UB,
                                     num_logs,
                                     lambdamat,
                                     keyvec,
                                     etavec,
                                     key0vec,
                                     eta0vec,
                                     alphavec,
                                     muvec,
                                     probvec_V,
                                     probvec_Z,
                                     probvec_Q,
                                     probvec_F,
                                     V_left,
                                     P0,
                                     minimum_n,
                                     open_segment,
                                     end_point,
                                     true);


  int V_12p = out_segment_12p["V"] ;  //V state for the new "merged" segment
  int Z_12p = out_segment_12p["Z"] ;  //Z state for the new "merged" segment
  int Q_12p = out_segment_12p["Q"] ;  //Q state for the new "merged" segment
  int F_12p = out_segment_12p["F"] ;  //F state for the new "merged" segment

  double sum_proposal_12p = out_segment_12p["sum_proposals"] ; // sum of log proposals for V_12p, Z_12p, Q_12p, F_12p
  double posterior_12p = out_segment_12p["sum_posterior"] ;  // log joint posterior for quantities involved in the new "left" segment



   // 4. evaluate states for the current "left" segment: V_1, Z_1, Q_1, F_1
   // 5. evaluate the proposal distribution at the current values
   // 6. evaluate the posterior of the current "left" segment


   List out_segment_1 = segment_box(Tvec,
                                     Yvec,
                                     V_1,
                                     Z_1,
                                     Q_1,
                                     F_1,
                                     U,
                                     K,
                                     W,
                                     LB,
                                     Td,
                                     num_logs,
                                     lambdamat,
                                     keyvec,
                                     etavec,
                                     key0vec,
                                     eta0vec,
                                     alphavec,
                                     muvec,
                                     probvec_V,
                                     probvec_Z,
                                     probvec_Q,
                                     probvec_F,
                                     V_left,
                                     P0,
                                     minimum_n,
                                     false,
                                     end_point,
                                     false);

  double sum_proposal_1 = out_segment_1["sum_proposals"] ; // sum of log proposals for V_1, Z_1, Q_1, F_1
  double posterior_1 = out_segment_1["sum_posterior"] ;  // log joint posterior for quantities involved in the current "left" segment


  // 7. evaluate states for the current "right" segment: V_2, Z_2, Q_2, F_2
  // 8. evaluate the proposal distribution at the current values
  // 9. evaluate the posterior of the current "right" segment


  List out_segment_2 = segment_box(Tvec,
                                   Yvec,
                                   V_2,
                                   Z_2,
                                   Q_2,
                                   F_2,
                                   U,
                                   K,
                                   W,
                                   Td,
                                   UB,
                                   num_logs,
                                   lambdamat,
                                   keyvec,
                                   etavec,
                                   key0vec,
                                   eta0vec,
                                   alphavec,
                                   muvec,
                                   probvec_V,
                                   probvec_Z,
                                   probvec_Q,
                                   probvec_F,
                                   V_1,
                                   P0,
                                   minimum_n,
                                   open_segment,
                                   end_point,
                                   false);

  double sum_proposal_2 = out_segment_2["sum_proposals"] ; // sum of log proposals for V_2, Z_2, Q_2, F_2
  double posterior_2 = out_segment_2["sum_posterior"] ;  // log joint posterior for quantities involved in the current "right" segment


  //logarithm of the Jacobian determinant l.jac

  double   l_jac = - std::log(L_12p)    ;

  // logarithm of the ratio btw the props of prosing a jump backward and  prosing a jump forward

  double log_jump = std::log(Jss1)-std::log(Js1s) ;

  // logarithm of the ratio involving the "right segment"

  double logr_rightV = 0 ;

  if(V_2 == 0 ){logr_rightV =  std::log(1-P0) ; }


  //final log ratio

    double log_ratio_b = posterior_12p - posterior_1 - posterior_2

                       + sum_proposal_1 + sum_proposal_2 - sum_proposal_12p

                       + l_jac + log_jump + logr_rightV ;


  //evaluate the move (toss a coin)

  double  coin = std::log( R::runif(0.0,1.0) ) ;


  if(coin<log_ratio_b){  //move is accepted: the interval is dropped

    // MOVE ACCEPTED: UPDATE THE STATE VECTORS Vvec, Zvec, Qvec, Fvec

    reduce_Bvec(Bvec, BR, S) ;

    reduce_vecZQFV(Zvec, Z_12p, BR, S) ;

    reduce_vecZQFV(Fvec, F_12p, BR, S) ;

    reduce_vecZQFV(Vvec, V_12p, BR, S) ;

    reduce_vecZQFV(Qvec, Q_12p, BR, S) ;

    S = S-1 ;


    NumericVector Bvec_actual = Bvec[Range(0,S)] ;
    bool check_order = std::is_sorted(Bvec_actual.begin(),Bvec_actual.end()) ;

    if(check_order == false){

      Rcout << "backward_function: this is Bvec:  " << Bvec_actual << "\n"  ;

      stop("backward_function: the Breakpoint vector is not sorted");

    }

    int illegal_Break = sum(Bvec_actual < start_point) ;

    if(illegal_Break > 1 ){

      Rcout << "start_point:  " << start_point << "\n" ;
      Rcout << "end_point:  " << end_point << "\n" ;
      Rcout << "backward_function: this is Bvec:  " << Bvec_actual << "\n"  ;
      stop("backward_function: more than 1 Breakpoint vector is smaller than LB");

      }

  }

  //end of the function

}





















