#include "Rcpp.h"
#include <algorithm>
#include "toolbox1.hpp"
#include "initial_state.hpp"
#include "initial_breakpoints.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' @title Function to initialize the MCMC in RJMCMC_SMC
//' Function generating the first iteration of the MCMC
//' @param T_seg Vector with time stamps within the update interval
//' @param Y_seg Vector with the messages within the update interval
//' @param minimum_n (int) minimum allowed number of observations in a active segments
//' @param start_point, (double)   the start point of the update interval
//' @param end_point, (double)   the end point of the update interval
//' @param t_star, (double) estimate of the last simulated breakpoint
//' @param K (int) number of Z states
//' @param W (int) Number of levels for the Q state
//' @param U (int) Number of V states
//' @param empty_mix (bool) If True (default) 2 kind of empty segments are assumed
//' @param probvec_V probability mass for the V state
//' @param probvec_Z probability mass for the Z state
//' @param probvec_Q probability mass for the Q state
//' @param probvec_F probability mass for the F state
//' @param alphavec Shape parameters ruling the rate of occurrence in different levels of Z state (non-empty segments)
//' @param muvec Mean parameters ruling the rate of occurrence in different levels of Z state
//' @param keyvec Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec Mean parameters ruling the length of empty segments in different levels of F state
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param lambdamat (const NumericMatrix&) each row v is the full prob vector for the messages of state V=v
//' @param num_logs Total number of different messages that can be observed (dictionary)
//' @param max_range (const double&) time interval in which the last changepoint is allowed to be place after the end point (end_point;end_point+1)
//' @param Smax maximum number of segments allowed inside the update interval
//' @return A list containting 8 elements:
//' \describe{
//' \item{S}{int, number of segments generarted}
//' \item{Bvec}{NumericVector, hosts the changepoints generated (inlcuding t_star, tau_k (last breakpoint after end_point)) }
//' \item{Vvec}{IntegerVector,  initial V state for each segment }
//' \item{Zvec}{IntegerVector,  initial Z state for each segment }
//' \item{Qvec}{IntegerVector,  initial Q state for each segment }
//' \item{Fvec}{IntegerVector,  initial F state for each segment }
//' \item{Mvec}{NumericVector,  initial M (rate) for each segment }
//' \item{Nvec}{NumericVector, count of observations in each segment }
//' }
//' @export
// [[Rcpp::export]]
List iteration0_RJMCMC (   const NumericVector& T_seg,
                           const IntegerVector& Y_seg,
                           const int& minimum_n,
                           const double& start_point,
                           const double& end_point,
                           const double& t_star,
                           const int& K,
                           const int& W,
                           const int& U,
                           const bool& empty_mix,
                           const NumericVector& probvec_V,
                           const NumericVector& probvec_Z,
                           const NumericVector& probvec_Q,
                           const NumericVector& probvec_F,
                           const NumericVector& alphavec,
                           const NumericVector& muvec,
                           const NumericVector& keyvec,
                           const NumericVector& etavec,
                           const NumericVector& key0vec,
                           const NumericVector& eta0vec,
                           const NumericMatrix lambdamat,
                           const double& P0,
                           const int& num_logs,
                           const double& max_range,
                           const int& Smax){


    // generate initial Breakpoint Vector

    double min_seg_length = 1/720 ;

    NumericVector Bvec = initial_breakpoints( T_seg,
                                              t_star,
                                              start_point,
                                              end_point,
                                              minimum_n,
                                              max_range,
                                              min_seg_length) ;


    // number of segments
    int S = Bvec.size() - 1 ;


    // Compute states for the generated segments


   List out_state= initial_state( T_seg,
                                  Y_seg,
                                  S,
                                  Bvec,
                                  end_point,
                                  probvec_V,
                                  probvec_Q,
                                  probvec_Z,
                                  probvec_F,
                                  lambdamat,
                                  U,
                                  W,
                                  K,
                                  alphavec,
                                  muvec,
                                  keyvec,
                                  etavec,
                                  key0vec,
                                  eta0vec,
                                  P0,
                                  num_logs,
                                  empty_mix,
                                  minimum_n);



     IntegerVector Vvec_partial =clone(as<IntegerVector>( out_state["Vvec"])) ;
     IntegerVector Zvec_partial =clone(as<IntegerVector>(  out_state["Zvec"]));
     IntegerVector Fvec_partial =clone(as<IntegerVector>(  out_state["Fvec"])) ;
     IntegerVector Qvec_partial =clone(as<IntegerVector>(  out_state["Qvec"])) ;
     IntegerVector Nvec_partial =clone( as<IntegerVector>( out_state["Nvec"])) ;


     if(S>Smax){ stop("number of segments S greater than Smax. Choose a bigger Smax");}

     NumericVector Bvec_final = NumericVector(Smax+1) ;
     Bvec_final[Range(0,S)] = Bvec ;

     IntegerVector Zvec_final = IntegerVector(Smax) ;
     Zvec_final[Range(0,S-1)] = Zvec_partial ;

     IntegerVector Vvec_final = IntegerVector(Smax) ;
     Vvec_final[Range(0,S-1)] = Vvec_partial ;

     IntegerVector Qvec_final = IntegerVector(Smax) ;
     Qvec_final[Range(0,S-1)] = Qvec_partial ;

     IntegerVector Fvec_final = IntegerVector(Smax) ;
     Fvec_final[Range(0,S-1)] = Fvec_partial ;



       return (List::create(
                           Named("Vvec")=Vvec_final,
                           Named("Zvec")=Zvec_final,
                           Named("Qvec")=Qvec_final,
                           Named("Fvec")=Fvec_final,
                           Named("S")=S,
                           Named("Nvec")=Nvec_partial,
                           Named("Bvec")=Bvec_final
                           ));

    }



























