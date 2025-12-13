#include <algorithm>
#include <RcppArmadillo.h>
#include "toolbox1.hpp"
#include "RJMCMC_SMC.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//' SMC Function implementing particle filter for the model in  Gramuglia et al. 2020 using the importance density approximation
//' as in Heard & Turcotte (2017)
//'
//' Given a time series of time-stamp and messages this function run a particle filter with
//' Reversible Jump Markov Chain Montecarlo steps for approximating the posterior distribution for the following variables
//' \itemize{
//'  \item the number of segments
//'  \item the position of the breakpoints defining the segments
//'  \item the state dynamic for each interval which include
//'  \itemize{
//'   \item Z, ruling the length of each non-empty segment
//'   \item Q, ruling the rate of occurrence of the messages in non-empty segments
//'   \item V, ruling the segment composition (in term of messages) in non-empty segments
//'   \item F, ruling the length of empty-segments
//'   }
//'}
//' The algorithm approximate these posterior distributions by splitting the observation
//' time interval in sub-intervals defined by the custom parameter \code{length_UI}
//' @param Yvec IntegerVector  with all the messages
//' @param Tvec NumericVector with all the time-stamps
//' @param length_UI length of the each update interval
//' @param n_particle number of particles created
//' @param U Number of levels of the V state
//' @param W Number of levels of the Q state
//' @param K Number of levels of the Z state
//' @param num_logs Total number of different messages that can be observed (dictionary)
//' @param lambdamat matrix U*num_logs with probability mass for the messages in each V state
//' @param keyvec Shape parameters ruling the length of non-empty segments in different levels of Z state
//' @param etavec Mean parameters ruling the length of non-empty segments in different levels of Z state
//' @param key0vec Shape parameters ruling the length of empty segments in different levels of F state
//' @param eta0vec Mean parameters ruling the length of empty segments in different levels of F state
//' @param alphavec Shape parameters ruling the rate of occurrence in different levels of Z state (non-empty segments)
//' @param muvec Mean parameters ruling the rate of occurrence in different levels of Z state
//' @param probvec_V probability mass for the V state
//' @param probvec_Z probability mass for the Z state
//' @param probvec_Q probability mass for the Q state
//' @param probvec_F probability mass for the F state
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param minimum_n minimum number of observations that must be observed in a non-empty segment
//' @param Jss1, (const double&) Probability to propose a jump forward
//' @param Js1s, (const double&) Probability to propose a jump backward
//' @param Smax maximum number of segments allowed inside the update interval
//' @param n_ite (const int&), number of iteration to be performed in the RJMCMC
//' @param burn_in (const int&) Number of iterations used as burn_in period.
//' @param thinning (const int&) 1 every "thinning" iterations will be used
//' @export
// [[Rcpp::export]]

  List SMC_turcotte_cpp( const IntegerVector& Yvec,
                         const NumericVector& Tvec,
                         const double& length_UI,
                         const int n_particle,
                         const int U,
                         const int W,
                         const int K,
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
                         const double& P0,
                         const int& minimum_n,
                         const double& Jss1,
                         const double& Js1s,
                         const int& Smax,
                         const int& n_ite,
                         const int& burn_in,
                         const int& thinning){


  //   @@@ START INITILIZATION STEP @@@

  double first_point =  floor( Tvec[0])   ;

  double last_point =  max(Tvec) + length_UI  ;

  NumericVector UI_bounds = as<NumericVector>(seq_cpp(first_point,last_point,length_UI) ); // boundaries for the Update Intervals

  NumericVector weight_vec = rep(1.0,n_particle);  // the normalized weight vector

  IntegerVector V_last_complete = rep(1,n_particle) ; // V state for the last complete segment of each particle

  double t_star_empty_seg = 0 ; //estimate of the last breakpoint at the end of each update interval for empty segments
  double t_star_non_empty_seg = 0 ; //estimate of the last breakpoint at the end of each update interval for non-empty segments

  // last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration
  NumericVector B_last = rep(0.0,n_particle) ;

  int n_UI = UI_bounds.size()-1 ; // number of Update intervals

  double start_point = UI_bounds[0] ; // left bound of the current Update Interval
  double end_point = 0.0 ; // right bound of the current Update Interval


  //initialize the containers holding the current Bvec, Vvec, Zvec, Qvec, Fvec, Svec at each SMC iteration

  List container_B(n_particle); // List for storing the Bvec generated inside the current Update Interval
  List container_V(n_particle); // List for storing the Vvec generated inside the current Update Interval
  List container_Z(n_particle); // List for storing the Zvec generated inside the current Update Interval
  List container_Q(n_particle); // List for storing the Qvec generated inside the current Update Interval
  List container_F(n_particle); // List for storing the Fvec generated inside the current Update Interval
  IntegerVector Svec(n_particle) ; // number of segments currently in each particle

  //initialize the final outputs

  List storage_B(n_UI); // List with the "container_B" generated at the end of each SMC iteration
  List storage_V(n_UI); // List with the "container_V" generated at the end of each SMC iteration
  List storage_Z(n_UI); // List with the "container_Z" generated at the end of each SMC iteration
  List storage_Q(n_UI); // List with the "container_Q" generated at the end of each SMC iteration
  List storage_F(n_UI); // List with the "container_F" generated at the end of each SMC iteration
  List storage_S(n_UI); // List with the "Svec" generated at the end of each SMC iteration
  List storage_weight(n_UI); // List with the normalize particle weights at the end of each SMC iteration



  // START SMC-RJMCMC


  for(int j=0; j<n_UI ; j++ ){

//Rcout << "@@@@@@@@@@@   updating window n  "  << j << "\n" ;

    end_point = UI_bounds[j+1] ;

    if(sum(Tvec>start_point & Tvec<end_point)==0){

      //Do nothing - no observations inside the current Update Interval

    }else{

      // select the observations to be used in the current Update Interval.
      // select all the obs recorded after min(B_last). This is because
      // the current posterior is based on the obs from "B_last" while the proposal from "t_star"
      // and being t_star ave(B_last) then min(B_last) < t_star is always true

      double B_last_min = min(B_last) ;
      double B_last_max = max(B_last) ;

      NumericVector  T_seg = Tvec[(Tvec>=B_last_min) & (Tvec<end_point)] ;
      IntegerVector  Y_seg = Yvec[(Tvec>=B_last_min) & (Tvec<end_point)] ;

      //perform RJMCMC

      // --- collect B_last value

      // retrieve largest_obs: largest observation between [B_last_min; B_last_max]

      NumericVector obs_in_B = T_seg[T_seg <= B_last_max] ;

      double largest_obs = max(obs_in_B) ;

      LogicalVector B_last_non_empty_seg_logical = B_last <=  largest_obs ;
      LogicalVector B_last_empty_seg_logical = B_last >  largest_obs ;


      NumericVector B_last_non_empty_seg = B_last[B_last_non_empty_seg_logical] ;
      NumericVector B_last_empty_seg =  B_last[B_last_empty_seg_logical] ;

      // compute "t_star for non-empty segments
      t_star_non_empty_seg = mean(B_last_non_empty_seg) ;

      // compute "t_star for empty segments
      t_star_empty_seg = mean(B_last_empty_seg) ;

      // indexes of the particles
      IntegerVector seq_particle = seq(0,n_particle-1) ;

      IntegerVector non_empty_particle_index_vec = seq_particle[B_last_non_empty_seg_logical] ;
      IntegerVector empty_particle_index_vec = seq_particle[B_last_empty_seg_logical] ;


      if(non_empty_particle_index_vec.length()>0){

        // RJMCMC SMC for non-empty segments
        RJMCMC_SMC(T_seg,
                   Y_seg,
                   U,
                   K,
                   W,
                   start_point,
                   end_point,
                   t_star_non_empty_seg,
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
                   P0,
                   minimum_n,
                   Jss1,
                   Js1s,
                   Smax,
                   n_ite,
                   burn_in,
                   thinning,
                   n_particle,
                   non_empty_particle_index_vec,
                   V_last_complete,
                   B_last,
                   container_B,
                   container_V,
                   container_Z,
                   container_Q,
                   container_F,
                   Svec,
                   weight_vec);
      }

      if(empty_particle_index_vec.length()>0){

        // RJMCMC SMC for empty segments
        RJMCMC_SMC(T_seg,
                   Y_seg,
                   U,
                   K,
                   W,
                   start_point,
                   end_point,
                   t_star_empty_seg,
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
                   P0,
                   minimum_n,
                   Jss1,
                   Js1s,
                   Smax,
                   n_ite,
                   burn_in,
                   thinning,
                   n_particle,
                   empty_particle_index_vec,
                   V_last_complete,
                   B_last,
                   container_B,
                   container_V,
                   container_Z,
                   container_Q,
                   container_F,
                   Svec,
                   weight_vec);
      }

      // normalize the weight vector

      double sum_weight = sum(weight_vec) ;

      for(int i=0; i<n_particle; i++){

        weight_vec[i] = weight_vec[i]  / sum_weight ;

      }

      //performe re-sampling (if needed)

      // // compute ESS value
      //
      // double ESS = 1.0 /sum( (Rcpp::pow(no_norm_weight_vec,2) )) ;
      //
      //
      // if(ESS < (M / 3) ){
      //
      //   /// @@@ START RESAMPLING THE PARTICLES GENERATED IN THIS ITERATIONS
      //
      //   resampling_window( weight_vec,
      //                      container_B,
      //                      container_Z,
      //                      container_V) ;
      //
      //   weight_vec = rep(1.0/M,M) ;
      //
      // }
      //store the weights
      //array_W[j] = clone(weight_vec) ;

      //store results
      //store the containers current Update Interval
      //
       storage_B[j] = clone( as<List>(container_B ) ) ;
       storage_V[j] = clone( as<List>(container_V ) ) ;
       storage_Z[j] = clone( as<List>(container_Z ) ) ;
       storage_Q[j] = clone( as<List>(container_Q ) ) ;
       storage_F[j] = clone( as<List>(container_F ) ) ;
       storage_S[j] = clone(Svec) ;
       storage_weight[j] = clone( as<List>(weight_vec ) ) ;

      start_point = end_point ;

    }




  }

  return ( List::create(Named("storage_B") = storage_B,
                        Named("storage_V") = storage_V,
                        Named("storage_Z") = storage_Z,
                        Named("storage_Q") = storage_Q,
                        Named("storage_F") = storage_F,
                        Named("storage_S") = storage_S,
                        Named("storage_weight") = storage_weight,
                        Named("n_UI") = n_UI,
                        Named("UI_bounds") = UI_bounds));




}










