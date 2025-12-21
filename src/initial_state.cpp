#include <algorithm>
#include <Rcpp.h>
#include "toolbox1.hpp"
#include "full_conditional_V.hpp"
#include "full_conditional_Z.hpp"
#include "collapsed_full_conditional_Q.hpp"
#include "full_conditional_F.hpp"
#include "extract_index.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


//' @title initial_state
//' Function for generating the initial state values in each segment
//' @param T_seg Vector with time stamps within the update interval
//' @param Y_seg Vector with the messages within the update interval
//' @param S Number of segments within the update interval
//' @param Bvec Position of the Changepoints simulated within the update interval
//' @param end_point (const double&) end point of the update interval
//' @param probvec_V probability mass for the V state
//' @param probvec_Q probability mass for the Q state
//' @param probvec_Z probability mass for the Z state
//' @param probvec_F probability mass for the F state
//' @param lambdamat matrix U*num_logs with probability mass for the messages in each V state
//' @param U Number of levels of the V state
//' @param W Number of levels of the Q state
//' @param K Number of levels of the Z state
//' @param alphavec Shape parameters ruling the rate of occurrence in different levels of Z state (non-empty segments)
//' @param muvec Mean parameters ruling the rate of occurrence in different levels of Z state
//' @param keyvec Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec Mean parameters ruling the length of empty segments in different levels of F state
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param num_logs Total number of different messages that can be observed (dictionary)
//' @param empty_mix Logical. Equal "true" if the empty segments are a mixture of the states
//' @param minimum_n minimum number of observations that must be observed in a non-empty segment
//' @return A list containting 6 elements:
//' \describe{
//' \item{Vvec}{IntegerVector,  initial V state for each segment }
//' \item{Zvec}{IntegerVector,  initial Z state for each segment }
//' \item{Qvec}{IntegerVector,  initial Q state for each segment }
//' \item{Fvec}{IntegerVector,  initial F state for each segment }
//' \item{Mvec}{NumericVector,  initial M (rate) for each segment }
//' \item{Nvec}{NumericVector, count of observations in each segment }
//' }
//' @export
// [[Rcpp::export]]
List initial_state ( const NumericVector& T_seg,
                     const IntegerVector& Y_seg,
                     const int& S,
                     const NumericVector& Bvec,
                     const double& end_point,
                     const NumericVector& probvec_V,
                     const NumericVector& probvec_Q,
                     const NumericVector& probvec_Z,
                     const NumericVector& probvec_F,
                     const NumericMatrix& lambdamat,
                     const int& U,
                     const int& W,
                     const int& K,
                     const NumericVector& alphavec,
                     const NumericVector& muvec,
                     const NumericVector& keyvec,
                     const NumericVector& etavec,
                     const NumericVector& key0vec,
                     const NumericVector& eta0vec,
                     const double& P0,
                     const int& num_logs,
                     const bool& empty_mix,
                     const int& minimum_n){



   // initialize the output vectors

    IntegerVector Vvec(S) ;      // V state vector
    IntegerVector Zvec(S) ;      // Z state vector
    IntegerVector Fvec(S) ;      // F state vector
    IntegerVector Qvec(S) ;      // Q state vector
     NumericVector Mvec(S) ;      // M rate vector



   // Lvec: vector with length of the segments
   NumericVector Lvec = diff(Bvec) ;

   //0) generate rate in each segment
   IntegerVector interval = cpp_findInterval(T_seg, Bvec) ;

   // Nvec: Vector with number of observations in each segment
   NumericVector Nvec = make_table(interval ,S,false) ;

   Mvec = Nvec / Lvec ;

   // V value for the previous segment (initialized at 1 so the first segment can be anything)
   int V_left = 1 ;


   //initialize open_segment: will be true when the last interval must be update

   bool open_segment = false;


  //sample V,Z,Q for the non-empty segments and F for the empty segments

  for(int i=0; i<S; i++){

     //int N_value = Nvec[i] ;


      if(Nvec[i]>0){   //Non-empty segment: sample V,Z,Q


         //if the last segment is selected, change the value of open_segment
         if(i==S-1){

            open_segment=true ;
         }

         //lower/upper bound for the specific segment
         double LB  = Bvec[i] ;
         double UB  = Bvec[i+1] ;

         // number of messages for the specific segment
         int N_seg = Nvec[i] ;

         // sample V
          List out_V = full_conditional_V(0,
                                          U,
                                          num_logs,
                                          T_seg,
                                          Y_seg,
                                          LB,
                                          UB,
                                          lambdamat,
                                          probvec_V,
                                          V_left,
                                          P0,
                                          end_point,
                                          open_segment,
                                          true);



         Vvec[i] = out_V["V"] ;

         // ENHANCED BUG DETECTION CHECK: Open segment (last segment) with zero observations must have V = 0
         // Directly count observations in the open segment [LB, end_point)
         int actual_count_open = 0;
         if(i == S-1){
           double LB_seg = Bvec[i];
           for(int k = 0; k < T_seg.size(); k++){
             if(T_seg[k] >= LB_seg && T_seg[k] < end_point){
               actual_count_open++;
             }
           }
           if((Nvec[i] == 0 || actual_count_open == 0) && Vvec[i] != 0){
             Rcpp::stop("BUG DETECTED in initial_state: Open segment (last segment) with zero observations has non-zero V state. Segment index=%d, V=%d, Nvec[%d]=%.0f, actual_count=%d, LB=%.6f, UB=%.6f, end_point=%.6f, T_seg_size=%d", 
                        i, Vvec[i], i, Nvec[i], actual_count_open, Bvec[i], Bvec[i+1], end_point, T_seg.size());
           }
         }

         V_left = out_V["V"] ;

         //sample Z
         double L_int =Lvec[i]  ;

         List out_Z = full_conditional_Z(0,
                                         K,
                                         L_int,
                                         keyvec,
                                         etavec,
                                         probvec_Z,
                                         true) ;


          Zvec[i] = as<int>(out_Z["Z"]) ;

         //sample Q

         List out_Q = collapsed_full_conditional_Q(0,
                                                   W,
                                                   N_seg,
                                                   LB,
                                                   UB,
                                                   alphavec,
                                                   muvec,
                                                   probvec_Q,
                                                   minimum_n,
                                                   open_segment,
                                                   end_point,
                                                   true);


          Qvec[i] = out_Q["Q"] ;

         // set F=0 for the current segment (non-empty segment)

           Fvec[i] = 0 ;


      }else{

         // sample F

         List out_F = full_conditional_F(0,
                                         Lvec[i],
                                         key0vec,
                                         eta0vec,
                                         probvec_F,
                                         true);


         Fvec[i] = out_F["F"] ;

         // set V=0, Z=0, Q=0 for the current segment (empty segment)

          Vvec[i] = 0 ;
          Zvec[i] = 0 ;
          Qvec[i] = 0 ;

          V_left = 0 ;

     }

  }




   // 5) return Vvec

   return( List::create(Named("Vvec")=Vvec,
                        Named("Nvec")=Nvec,
                        Named("Zvec")=Zvec,
                        Named("Mvec")=Mvec,
                        Named("Qvec")=Qvec,
                        Named("Fvec")=Fvec
                           ));



 }





