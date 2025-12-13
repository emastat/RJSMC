#include <algorithm>
#include <Rcpp.h>
#include "segment_box.hpp"
#include "extract_index.hpp"
#include "updating_vec.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @title shift_function: This function takes in the current state of the MCMC and employ a position shift of a breakpoint.
//' all the input parameters described below refers to their most recent values.
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

// [[Rcpp::export]]

void  shift_function(      int& S,
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

    //Rcout << "start shift " << "\n" ;

    NumericVector Bvec_clone = clone(Bvec) ;


    int BR ; // index (from 1 to S-1) of the breakpoint to shifted

    double T_p=0, // the proposed changepoint
           u1=0, // value of a Unif(0,1) used to sampled a new changepoint
           LB=0,  // lower bound of the selected "double" segment
           UB=0,  // upper bound of the selected "double" segment
           Td=0,  // the breakpoint that will be possibly shifted
           L_1p=0, //length of the new "left" segment
           L_2p=0; //length of the new "right" segment

    int V_1 = 0, // V state of the current "left" segment
        V_2 = 0, // V state of the current "right" segment
        Z_1 = 0, // Z state of the current "left" segment
        Z_2 = 0, // Z state of the current "right" segment
        Q_1 = 0, // Q state of the current "left" segment
        Q_2 = 0, // Q state of the current "right" segment
        F_1 = 0, // F state of the current "left" segment
        F_2 = 0, // F state of the current "right" segment
        V_left = 0, //V  state of the segment before the current one
        V_right = 0 ; //V state of the segment after the current one

    double count_1p = 0.0 ; // number of obs in proposed seg 1
    double count_2p = 0.0 ; // number of obs in proposed seg 2

    int counter = 0 ;

    bool open_segment = false ; // to check if the last segment has been proposed to be changed

    bool illegal_configuration = true ; // to check that the split provides a legal confguration of the new segments (given te assumed constrains)


    while( (illegal_configuration==true) & (counter<50)){

      L_1p=0,
      L_2p=0;

      while((L_1p< (1.0/720.0)) | (L_2p< (1.0/720.0))){

        BR =as<int>( Rcpp::sample(S-1, 1)) ;

        Td = Bvec[BR] ;

        LB = Bvec[BR-1]  ;

        UB = Bvec[BR+1]  ;

        u1 =  R::runif(0.0,1.0) ;

        if((BR>1) & (BR < (S-1))){ // the new changepoint can be placed anywhere within the selected segment

           T_p =  LB + u1*(UB-LB)  ;

        }else if((BR>1) & (BR==(S-1))){ // the new changepoint must be placed within the penultimate changepoint (LB) and the end of the Update interval (end_point)

           T_p = LB + u1*(end_point-LB)  ;

        }else if((BR==1) & (BR<(S-1))){ // the new changepoint must be placed within the start point of the update interval start_point) and the 3rd changepoint (UB)

           T_p = start_point + u1 * (UB - start_point) ;

        }else if((BR==1) & (BR==(S-1))){ //the new changepoint must be placed between start_point and end_point

           T_p = start_point + u1 * (end_point - start_point) ;
        }

        L_1p = T_p - LB ;

        L_2p = UB - T_p ;

        counter += 1 ;

        if(counter>50){ break ;}
      }

      // if the last segment has been sample then it must be an open segment

      if(BR == (S-1)){
        open_segment = true ;

        if(T_p>end_point){

          Rcout << "LB  " << LB << "\n" ;
          Rcout << "T_p  " << T_p << "\n" ;
          Rcout << "UB  " << UB << "\n" ;
          Rcout << "end_point  " << end_point << "\n" ;
          stop("the last segment must end after the end_point") ;
        }

      }else{ open_segment = false ;}

      // to compute the move it is required;

      V_1 = Vvec[BR-1]  ;
      V_2 = Vvec[BR] ;

      Z_1 = Zvec[BR-1]  ;
      Z_2 = Zvec[BR]   ;

      Q_1 = Qvec[BR-1]  ;
      Q_2 = Qvec[BR]   ;

      F_1 = Fvec[BR-1]  ;
      F_2 = Fvec[BR]   ;

      if( BR>1 ){   V_left = Vvec[BR-2] ; }else{ V_left = -1 ; }


      if( BR<(S-1)){    V_right = Vvec[BR+1] ; }else{  V_right = -1 ; }


      V_1 = Vvec[BR-1]  ;

      V_2 = Vvec[BR]  ;

      Z_1 = Zvec[BR-1]  ;

      Z_2 = Zvec[BR]  ;

      Q_1 = Qvec[BR-1]  ;

      Q_2 = Qvec[BR]  ;

      F_1 = Fvec[BR-1]  ;

      F_2 = Fvec[BR]  ;


      NumericVector out_extract_1p =clone( extract_index(Tvec,LB,T_p)) ;

      NumericVector out_extract_2p(3) ;

      if(BR < (S-1) ){

         out_extract_2p =clone( extract_index(Tvec,T_p,UB)) ;

      }else{

         out_extract_2p =clone( extract_index(Tvec,T_p,end_point)) ;

      }


      count_1p = out_extract_1p[0] ; // number of obs in [LB;T_p), the new "left" segment
      count_2p = out_extract_2p[0] ; // number of obs in [T_p;UB), the new "rigt" segment


      if(( (count_1p==0) & (V_left == 0) ) |
         ( (count_2p==0) & (V_right == 0) ) |
         ( (count_1p>0) & (count_1p<minimum_n) ) |
         ( (count_2p>0) & ( count_2p<minimum_n) )  ){


         //the new proposed changepoint produces two empty segments in a row
         //Illegal configuration stays "true": propose a new changepoint

      }else{

        // a legal configuration has been proposed, we can move on and performing the shift move

        illegal_configuration =false ;

      }

    }


   if(counter>49){

      // exit the forward function; too small segments are created

   }else{

   // 1. sample states for the  new "left" segment: V_1p, Z_1p, Q_1p, F_1p
   // 2. evaluate the proposal distribution at the new sampled values
   // 3. evaluate the posterior of the new "left" segment

   if(( (count_1p > 0) & (count_1p < minimum_n) ) |
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
      stop("error in shift: the 2 proposed segments have n_obs within 0 and mininum_n and this cannot be ");

   }

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


   // 7. evaluate states for the current "left" segment: V_1, Z_1, Q_1, F_1
   // 8. evaluate the proposal distribution at the current values
   // 9. evaluate the posterior of the current "left" segment

   List out_segment_1  = segment_box(Tvec,
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
   double posterior_1 = out_segment_1["sum_posterior"] ;  // log joint posterior for quantities involved in the current  "left" segment


   // 10. evaluate states for the current "right" segment: V_2, Z_2, Q_2, F_2
   // 11. evaluate the proposal distribution at the current values
   // 12. evaluate the posterior of the current "right" segment

   List out_segment_2  = segment_box(Tvec,
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
                                     V_left,
                                     P0,
                                     minimum_n,
                                     open_segment,
                                     end_point,
                                     false);


   double sum_proposal_2 = out_segment_2["sum_proposals"] ; // sum of log proposals for V_2, Z_2, Q_2, F_2
   double posterior_2 = out_segment_2["sum_posterior"] ;  // log joint posterior for quantities involved in the current "right" segment



     // logarithm of the ratio involving the "right segment"

     double logr_rightV = 0 ;

     if((V_2!= 0) & (V_2p==0)  ){logr_rightV = - std::log(1-P0) ; }

     else if( (V_2 ==0) & (V_2p !=0)){logr_rightV = + std::log(1-P0) ; }

   //final log ratio


   double log_ratio_s = posterior_1p + posterior_2p - posterior_1 - posterior_2

                      + sum_proposal_1 + sum_proposal_2 - sum_proposal_1p - sum_proposal_2p

                      + logr_rightV ;


   //evaluate the move (toss a coin)

   double  coin = std::log( R::runif(0,1) ) ;


    if( coin < log_ratio_s){

       // MOVE ACCEPTED: UPDATE THE STATE VECTORS Vvec, Zvec, Qvec, Fvec


       Bvec[BR] = T_p ;

       Vvec[BR-1] = V_1p ;
       Vvec[BR]   = V_2p ;

       Zvec[BR-1] = Z_1p ;
       Zvec[BR]   = Z_2p ;

       Qvec[BR-1] = Q_1p ;
       Qvec[BR]   = Q_2p ;

       Fvec[BR-1] = F_1p ;
       Fvec[BR]   = F_2p ;


       NumericVector Bvec_actual = Bvec[Range(0,S)] ;
       bool check_order = std::is_sorted(Bvec_actual.begin(),Bvec_actual.end()) ;

       if(check_order == false){

          Rcout << "shift_function: this is Bvec:  " << Bvec_actual << "\n"  ;

          stop("shift_function: the Breakpoint vector is not sorted");

       }

       int illegal_Break = sum(Bvec_actual < start_point) ;

       if(illegal_Break > 1 ){

          Rcout << "start_point:  " << start_point << "\n" ;
          Rcout << "end_point:  " << end_point << "\n" ;
          Rcout << "original_bvec:  " << Bvec_clone << "\n" ;
          Rcout << "shift_function: this is Bvec:  " << Bvec_actual << "\n"  ;
          stop("shit_function: more than 1 Breakpoint vector is smaller than LB");

       }


    }

   }
    //end of the function

 }













