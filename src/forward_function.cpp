#include <algorithm>
#include <Rcpp.h>
#include "segment_box.hpp"
#include "extract_index.hpp"
#include "updating_vec.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @title forward_function
//' This function takes in the current state of the MCMC and employ a forward step where a split of an interval is proposed
//' all the input parameters described below refers to their  most recent values.
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

void         forward_function(            int& S,
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


  //Rcout << "start forward " << "\n" ;

  //initialize quantities needed for performing the forward move

  int IN = 0 ;   // index (from 1 to S) of the segment to split

  double LB=0.0, // lower bound of the selected segment
         UB=0.0, // upper bound of the selected segment
         u1=0.0, // value of a Unif(0,1) used to sampled a new changepoint
         T_p=0.0, // the proposed changepoint
         L_12 = 0.0, // length of the current segment
         L_1p=0.0,  //length of the new  left segment
         L_2p=0.0;  // length of the second right segment

  int V_12 = 0, // V state of the current segment
      Z_12 = 0, // Z state of the current segment
      Q_12 = 0, // Q state of the current segment
      F_12 = 0, // F state of the current segment
      V_left = 0, //V  state of the segment before the current one
      V_right = 0 ; //V state of the segment after the current one

  double count_1p = 0.0 ; // number of obs in proposed seg 1
  double count_2p = 0.0 ; // number of obs in proposed seg 2

  bool open_segment = false ; // to check if the last segment has been proposed to be split

  int counter = 0 ; // count how many times an interval smaller than 0.5 second is met

  bool illegal_configuration = true ; // to check that the split provides a legal configuration of the new segment (given the assumed constrains)

  int count_illegal_configurations = 0 ; // count how many illegal configurations have been proposed

  bool too_many_small_seg = false; //flag equal to True if too many times (50) small segment has been proposed

  while(illegal_configuration==true){

    L_1p=0.0;  //length of the new  left segment
    L_2p=0.0;  // length of the second right segment


    //generate IN:  (the new segments must be bigger than 5 sec)

    while((L_1p < (1.0/720.0)) | (L_2p < (1.0/720.0))){


      //sample an interval to split

      IN =as<int>( Rcpp::sample(S, 1)) ;

      //lower and upper bound of the selected interval

      LB =  Bvec[IN-1]  ;

      UB =  Bvec[IN]  ;

      u1 =  R::runif(0.0,1.0) ;

      // generate the proposal breakpoint T_p

      if((IN>1) & (IN < S)){

        // the new changepoint can be placed anywhere within the selected segment
        // BUT: if IN > 0, T_p must be >= start_point (only Bvec[0] can be < start_point)
        // So we clamp LB to be at least start_point

        double effective_LB = std::max(LB, start_point);
        T_p =  effective_LB + u1*(UB-effective_LB)  ;


      }else if((IN==S) & (S>1)){

        // the new changepoint must be placed within the penultimate changepoint (LB) and the end of the Update interval (end_point)
        // BUT: if IN > 0, T_p must be >= start_point
        double effective_LB = std::max(LB, start_point);
        T_p =  effective_LB + u1*(end_point-effective_LB)  ;

      }else if((IN==1) & (S==1)){

        // there is only segment, so we can split it wherever withing the update interval
        // T_p must be >= start_point (IN==1 means we're inserting at position 1, not 0)

        T_p = start_point + u1*(end_point - start_point) ;

      }else if((IN==1) & (S>1)){

        // the new changepoint must be placed within the start point of the update interval (start_point) and the 3rd changepoint (UB)
        // T_p must be >= start_point (IN==1 means we're inserting at position 1, not 0)

        T_p = start_point + u1 * (UB - start_point) ;

      }

      L_1p = T_p - LB ;

      L_2p = UB - T_p ;

      counter += 1 ;

      if(counter>50){

        // too many small segment have been proposed, break the the loop

        too_many_small_seg = true ;
        break ;

      }

    }

    if(too_many_small_seg == true){

      // too many small segments, exit the forward function
      break;

    }else{

      // reset the counter - the proposed new interval are big enough
      counter = 0 ;

      // if the last segment has been sampled then it must be an open segment

      if(IN == S){

        open_segment = true ;

        if((T_p < LB) | (T_p > end_point)){

          Rcout << "LB:  " << LB << "\n" ;
          Rcout << "T_p:  " << T_p << "\n" ;
          Rcout << "UB:  " << UB << "\n" ;
          Rcout << "end_point:  " << end_point << "\n" ;
          Rcout << "IN: " << IN << "\n" ;
          Rcout << "S: " << S << "\n" ;
          Rcout << "u1:  " << u1 << "\n" ;

          stop("T_p must be btw LB & end_point") ;

        }

      }else{ open_segment = false ;}

      // to compute the move it is required;

      V_12 = Vvec[IN-1] ;

      Z_12 = Zvec[IN-1] ;

      Q_12 = Qvec[IN-1] ;

      F_12 = Fvec[IN-1];

      if( IN>1 ){   V_left = Vvec[IN-2] ; }else{ V_left = -1 ; }

      if( IN<S ){    V_right = Vvec[IN] ; }else{  V_right = -1 ; }

      L_12 = UB-LB ;

      if(L_12 == 0.0){stop("a L=0 segment is sampled in forward") ;}


      NumericVector out_extract_1p =clone( extract_index(Tvec,LB,T_p)) ;

      NumericVector out_extract_2p(3) ;
      NumericVector out_extract_12(3) ;

      if(IN < S){

        out_extract_2p =clone( extract_index(Tvec,T_p,UB)) ;

        out_extract_12 =clone( extract_index(Tvec,LB,UB)) ;

      }else{

        out_extract_2p =clone( extract_index(Tvec,T_p,end_point)) ;

        out_extract_12 =clone( extract_index(Tvec,LB,end_point)) ;

     }


      // number of obs in [LB;T_p), the new left segment
      count_1p = out_extract_1p[0] ;

      // number of obs in [T_p;UB), the new right segment
      count_2p = out_extract_2p[0] ;

      bool condition_1 = V_12==0 ;
      bool condition_2 = (count_1p==0) & (V_left == 0) ;
      bool condition_3 = (count_2p==0) & (V_right == 0) ;
      bool condition_4 = (count_1p>0)  & (count_1p<minimum_n) ;
      bool condition_5 = (count_2p>0) & (count_2p<minimum_n) & (open_segment==false) ;

      if( ( condition_1 ) |
          ( condition_2 ) |
          ( condition_3 ) |
          ( condition_4 ) |
          ( condition_5 )
        ) {

        //Illegal configuration stays "true": propose a new split
        count_illegal_configurations +=1 ;

      }else{

        // a legal configuration has been proposed, we can move on and performing the forward move

        illegal_configuration = false ;
      }

    }

    if(count_illegal_configurations > 49){

      // break the while loop, too many illegal configurations
      break ;

    }

  }

  if((count_illegal_configurations> 49) | (too_many_small_seg ==true)){

    // exit the forward function, not legal proposals have been generdted



  }else{


    if( ((count_1p > 0) & (count_1p < minimum_n)) |
       (( (count_2p > 0) & (count_2p < minimum_n)) & (open_segment == false)))
    {

      Rcout << "count_1p: " << count_1p << "\n" ;
      Rcout << "count_2p: " << count_2p << "\n" ;
      Rcout << "Tvev:" << Tvec << "\n" ;
      Rcout << "LB:" << LB << "\n" ;
      Rcout << "T_p:" << T_p << "\n" ;
      Rcout << "UB:" << UB << "\n" ;
      Rcout << "end_point:" << end_point << "\n" ;
      Rcout << "open_segment: " << open_segment  << "\n" ;
      stop("error in forward: the 2 proposed segments have n_obs within 0 and mininum_n and this cannot be ");

    }

    // 1. sample states for the  new "left" segment: V_1p, Z_1p, Q_1p, F_1p
    // 2. evaluate the proposal distribution at the new sampled values
    // 3. evaluate the posterior of the new "left" segment

    List out_segment_1p = segment_box(Tvec,
                                      Yvec,
                                      0,
                                      0,
                                      0,
                                      0,
                                      U,
                                      K,
                                      W,
                                      LB,
                                      T_p,
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
                                      true);


    int V_1p = out_segment_1p["V"] ;  //V state for the new "left" segment
    int Z_1p = out_segment_1p["Z"] ;  //Z state for the new "left" segment
    int Q_1p = out_segment_1p["Q"] ;  //Q state for the new "left" segment
    int F_1p = out_segment_1p["F"] ;  //F state for the new "left" segment

    double sum_proposal_1p = out_segment_1p["sum_proposals"] ; // sum of log proposals for V_1p, Z_1p, Q_1p, F_1p
    double posterior_1p = out_segment_1p["sum_posterior"] ;  // log joint posterior for quantities involved in the new "left" segment


    // 4. sample states for the new "right" segment: V_2p, Z_2p , Q_2p, F_2p
    // 5. evaluate the proposal distribution at the new sampled values
    // 6. evaluate the posterior of the new "right" segment

    List out_segment_2p = segment_box(Tvec,
                                      Yvec,
                                      0,
                                      0,
                                      0,
                                      0,
                                      U,
                                      K,
                                      W,
                                      T_p,
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
                                      V_1p,
                                      P0,
                                      minimum_n,
                                      open_segment,
                                      end_point,
                                      true);


    int V_2p = out_segment_2p["V"] ;  //V state for the new "right" segment
    int Z_2p = out_segment_2p["Z"] ;  //Z state for the new "right" segment
    int Q_2p = out_segment_2p["Q"] ;  //Q state for the new "right" segment
    int F_2p = out_segment_2p["F"] ;  //F state for the new "right" segment

    double sum_proposal_2p = out_segment_2p["sum_proposals"] ; // sum of log proposals for V_2p, Z_2p, Q_2p, F_2p
    double posterior_2p = out_segment_2p["sum_posterior"] ;  // log joint posterior for quantities involved in the new "right" segment

    // 7. evaluate states for the current segment: V_12, Z_12, Q_12, F_12
    // 8. evaluate the proposal distribution at the current values
    // 9. evaluate the posterior of the current segment

    List out_segment_12 = segment_box(Tvec,
                                      Yvec,
                                      V_12,
                                      Z_12,
                                      Q_12,
                                      F_12,
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
                                      false);


    double sum_proposal_12 = out_segment_12["sum_proposals"] ; // sum of log proposals for V_12, Z_12, Q_12, F_12
    double posterior_12 = out_segment_12["sum_posterior"] ;  // log joint posterior for quantities involved in the current segment

    //l_jac_ logarithm of the Jacobian determinant

    double l_jac = std::log(std::fabs(L_12)) ;

    //logarithm of the ratio btw the props of prosing a jump backward and  prosing a jump forward

    double log_jump = std::log(Js1s) - std::log(Jss1) ;


    // logarithm of the ratio involving the new "right" segment

    double logr_rightV = 0 ;

    if((V_2p==0) & (V_right!=0)){logr_rightV = - std::log(1-P0) ; }

    //final log ratio

    double log_ratio_f = posterior_1p + posterior_2p - posterior_12

                       + sum_proposal_12 - sum_proposal_1p - sum_proposal_2p

                       + l_jac + log_jump + logr_rightV   ;


    //evaluate the move (toss a coin)

    double  coin = std::log( R::runif(0,1) ) ;



    if( coin< log_ratio_f){

      // MOVE ACCEPTED: UPDATE THE STATE VECTORS Vvec, Zvec, Qvec, Fvec

      augment_Bvec(Bvec, T_p, IN, S ) ;

      augment_vecZQFV(Zvec, Z_1p, Z_2p, IN, S) ;

      augment_vecZQFV(Fvec, F_1p, F_2p, IN, S) ;

      augment_vecZQFV(Vvec, V_1p, V_2p, IN, S) ;

      augment_vecZQFV(Qvec, Q_1p, Q_2p, IN, S) ;

      S = S +1 ;

      NumericVector Bvec_actual = Bvec[Range(0,S)] ;
      bool check_order = std::is_sorted(Bvec_actual.begin(),Bvec_actual.end()) ;

      if(check_order == false){

        Rcout << "forward_function: this is Bvec:  " << Bvec_actual << "\n"  ;

        stop("forward_function: the Breapoint vector is not sorted");

      }

      int illegal_Break = sum(Bvec_actual < start_point) ;

      if(illegal_Break > 1 ){
        stop("forward_function: more than 1 Breakpoint vector is smaller than LB");
      }

    }


  }
//end of the function
}



















