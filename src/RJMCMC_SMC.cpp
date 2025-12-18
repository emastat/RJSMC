#include <Rcpp.h>
#include "iteration0_RJMCMC.hpp"
#include "forward_function.hpp"
#include "backward_function.hpp"
#include "shift_function.hpp"
#include "incremental_weight.hpp"
#include "weigths_correction.hpp"
#include "extract_index.hpp"
#include "toolbox2.hpp"

using namespace Rcpp;


//' @title RJMCMC_SMC:  function for generating samples from the important density through RJMCMC updates
//' @param T_seg Vector with time stamps within the update interval
//' @param Y_seg Vector with the messages within the update interval
//' @param U Number of levels of the V state
//' @param W Number of levels of the Q state
//' @param K Number of levels of the Z state
//' @param start_point start point of the update interval
//' @param end_point (const double&) end point of the update interval
//' @param t_star Estimate of the latest changepoint simulated when performing the RJMCMC during the previous iteration of the SMC (precedent update interval)
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
//' @param n_particle number of particles created
//' @param particle_index_vec vector with the index of the particles that will be updated
//' @param V_last_complete Vector with the V state for the last complete segment of each particle
//' @param B_last (NumericVector&) Vector with the last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration
//' @param container_B List for with Bvec generated inside the current Update Interval, for each particle
//' @param container_V List for with Vvec generated inside the current Update Interval, for each particle
//' @param container_Z  List for with Zvec generated inside the current Update Interval, for each particle
//' @param container_Q  List for with Qvec generated inside the current Update Interval, for each particle
//' @param container_F  List for with Fvec generated inside the current Update Interval, for each particle
//' @param Svec vector with number of segment in each particle
//' @param weight_vec vector with the current normalized weight for each particle
//' @param empty_seg (Logical) if true then the importance density for empty segments is used
//' @return The function does not explicitly return any object: However, by reference 9 element are updated:
//' \describe{
//' \item{container_B}{List for with Bvec generated inside the current Update Interval, for each particle}
//' \item{container_V}{List for with Vvec generated inside the current Update Interval, for each particle}
//' \item{container_Z}{List for with Zvec generated inside the current Update Interval, for each particle}
//' \item{container_Q}{List for with Qvec generated inside the current Update Interval, for each particle}
//' \item{container_F}{List for with Fvec generated inside the current Update Interval, for each particle}
//' \item{Svec}{vector with number of segment in each particle}
//' \item{B_last}{Vector with the last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration}
//' \item{V_last_complete}{Vector with the V state for the last complete segment of each particle}
//' \item{weight_vec}{vector with the updated  non normalized weight for each particle}
//' }
//' @export
// [[Rcpp::export]]
void     RJMCMC_SMC(const NumericVector& T_seg,
                    const IntegerVector& Y_seg,
                    const int& U,
                    const int& K,
                    const int& W,
                    const double& start_point,
                    const double& end_point,
                    const double& t_star,
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
                    const int& thinning,
                    const int& n_particle,
                    const IntegerVector& particle_index_vec,
                          IntegerVector& V_last_complete,
                          NumericVector& B_last,
                          List& container_B,
                          List& container_V,
                          List& container_Z,
                          List& container_Q,
                          List& container_F,
                          IntegerVector& Svec,
                          NumericVector& weight_vec,
                    const bool& empty_seg){


  //open_segment: true if in a RJMCMC iteration only 1 segment is generated (S=1), false otherwise

  bool open_segment = false ;

  IntegerVector particle_index_vec_local = clone((particle_index_vec));

  //randomly permute the indexes (permutation as in Heard & Turcotte, 2017)
  std::random_shuffle(particle_index_vec_local.begin(),particle_index_vec_local.end()) ;


  int particle_count = particle_index_vec_local.length() ;  //number of particles to update

  int particle_index = 0 ; //initialize index for particle currently updated

  int index = 0 ;      //initialize general index

  int V_first = 0 ; //V state for the first segment of the new simulated vector (for particle "i")


// 1. initialize the RJMCMC performing the initial iteration


  List out_iteration_0 = iteration0_RJMCMC( T_seg,
                                            Y_seg,
                                            minimum_n,
                                            start_point,
                                            end_point,
                                            t_star,
                                            K,
                                            W,
                                            U,
                                            true,
                                            probvec_V,
                                            probvec_Z,
                                            probvec_Q,
                                            probvec_F,
                                            alphavec,
                                            muvec,
                                            keyvec,
                                            etavec,
                                            key0vec,
                                            eta0vec,
                                            lambdamat,
                                            P0,
                                            num_logs,
                                            0.2,
                                            Smax);

    // check the all element of B_last are less than start_point
    for(int i=0; i<particle_index_vec_local.length(); i++){
      int ii = particle_index_vec_local[i];
      if(B_last[ii] > start_point){

        Rcout << "this is start_point   " << start_point << "\n" ;
        Rcout << "this is ii   " << ii << "\n" ;
        Rcout << "this is B_last[ii]   " << B_last[ii] << "\n" ;
        Rcpp::stop("RJMCMC_SMC.cpp at right after iteration 0: B_last[%d] = %f must be <= start_point = %f", 
                   ii, B_last[ii], start_point);
      }
    }                                            


    int S =  out_iteration_0["S"];

    // Bvec contains vector of breakpoints in the current Update Interval including t_star and tau_k (last breakpoint after end_point)
    NumericVector Bvec = out_iteration_0["Bvec"];

    IntegerVector Vvec = out_iteration_0["Vvec"];

    IntegerVector Zvec = out_iteration_0["Zvec"];

    IntegerVector Qvec = out_iteration_0["Qvec"];

    IntegerVector Fvec = out_iteration_0["Fvec"];

    int S_s ;

    IntegerVector urn = seq(1,3)  ;

    NumericVector jump_probs =NumericVector::create(Js1s,1-Js1s-Jss1,Jss1) ;


    if(n_ite<= burn_in){stop("n_ite cannot be smaller than burn_in");}

    // test for Bvec[1] < start_point                                        
    
    if( (S>1) & (Bvec[1] < start_point)){
      Rcpp::stop("RJMCMC_SMC.cpp at iteration 0: Bvec[1] = %f must be >= start_point = %f", 
                 Bvec[1], start_point);
    }

    if (Bvec[0] > start_point){
      Rcpp::stop("RJMCMC_SMC.cpp at iteration 0: Bvec[0] = %f must be <= start_point = %f", 
                 Bvec[0], start_point);
    }


     //2. perform  RJMCMC

    for(int i=0 ;i<n_ite ; i++){


      open_segment = false ;

      if((S>=2) & (S<Smax)){  S_s = as<int>( Rcpp::sample(urn, 1,false,jump_probs ) ) ;

      }else if(S==Smax ){  S_s = 1 ;

      }else if(S==1 ){S_s = 3 ;}


      if(S_s==3){

         forward_function(S,
                          T_seg,
                          Y_seg,
                          U,
                          K,
                          W,
                          probvec_V,
                          probvec_Z,
                          probvec_Q,
                          probvec_F,
                          lambdamat,
                          keyvec,
                          etavec,
                          key0vec,
                          eta0vec,
                          alphavec,
                          muvec,
                          num_logs,
                          Bvec,
                          Zvec,
                          Qvec,
                          Vvec,
                          Fvec,
                          Jss1,
                          Js1s,
                          P0,
                          minimum_n,
                          start_point,
                          end_point);

      // test for Bvec[0] < start_point                                        
      if( (S>1) & (Bvec[1] < start_point)){
        Rcpp::stop("RJMCMC_SMC.cpp at forward_function: Bvec[1] = %f must be >= start_point = %f", 
                   Bvec[1], start_point);
      }
      if (Bvec[0] > start_point){
        Rcpp::stop("RJMCMC_SMC.cpp at forward_function: Bvec[0] = %f must be <= start_point = %f", 
                   Bvec[0], start_point);
      }

      }else if(S_s==1){

         backward_function(S,
                          T_seg,
                          Y_seg,
                          U,
                          K,
                          W,
                          probvec_V,
                          probvec_Z,
                          probvec_Q,
                          probvec_F,
                          lambdamat,
                          keyvec,
                          etavec,
                          key0vec,
                          eta0vec,
                          alphavec,
                          muvec,
                          num_logs,
                          Bvec,
                          Zvec,
                          Qvec,
                          Vvec,
                          Fvec,
                          Jss1,
                          Js1s,
                          P0,
                          minimum_n,
                          start_point,
                          end_point);


      // test for Bvec[0] < start_point                                        
      if( (S>1) & (Bvec[1] < start_point)){
        Rcpp::stop("RJMCMC_SMC.cpp at backward_function: Bvec[1] = %f must be >= start_point = %f", 
                   Bvec[1], start_point);
      }

      if (Bvec[0] > start_point){
        Rcpp::stop("RJMCMC_SMC.cpp at backward_function: Bvec[0] = %f must be <= start_point = %f", 
                   Bvec[0], start_point);
      }


      }else{


            shift_function(S,
                           T_seg,
                           Y_seg,
                           U,
                           K,
                           W,
                           probvec_V,
                           probvec_Z,
                           probvec_Q,
                           probvec_F,
                           lambdamat,
                           keyvec,
                           etavec,
                           key0vec,
                           eta0vec,
                           alphavec,
                           muvec,
                           num_logs,
                           Bvec,
                           Zvec,
                           Qvec,
                           Vvec,
                           Fvec,
                           Jss1,
                           Js1s,
                           P0,
                           minimum_n,
                           start_point,
                           end_point);

      // test for Bvec[0] < start_point                                        
      if( (S>1) & (Bvec[1] < start_point)){
        Rcpp::stop("RJMCMC_SMC.cpp at shift_function: Bvec[1] = %f must be >= start_point = %f", 
                   Bvec[1], start_point);
      }

      if (Bvec[0] > start_point){
        Rcpp::stop("RJMCMC_SMC.cpp at shift_function: Bvec[0] = %f must be <= start_point = %f", 
                   Bvec[0], start_point);
      }
  
      } //end of the S_s condition

      V_first =Vvec[0] ;

      //if we are at burn in do nothing. Otherwise update and store the new state vectors

      if(i<burn_in){

      }else{

        //conditions:
        // 1. consider only 1 every "thinning" iterations
        // 2. there are still particles to update

        if( (i % thinning ==0 ) & (particle_count>0)  ){

          //select the index of the particles to update

          index = particle_index_vec_local[particle_index] ;

          // if all constraints are satisfied update the particle


          if(S==1){open_segment = true ; }

          // out_weight: the log incremental weight for the selected particle

          double out_weight = incremental_weight(T_seg,
                                                 Y_seg,
                                                 Vvec[0],
                                                 Zvec[0],
                                                 Qvec[0],
                                                 Fvec[0],
                                                 U,
                                                 K,
                                                 W,
                                                 B_last[index],
                                                 Bvec[1],
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
                                                 V_last_complete[index],
                                                 P0,
                                                 minimum_n,
                                                 t_star,
                                                 end_point,
                                                 open_segment);

          // START CODE TO REMOVE
          if(NumericVector::is_na(weight_vec[index])){

            stop("RJMCMC_SMC.cpp. weight_vec[index] is na. probably from before");

          }

          if(Rcpp::traits::is_nan<REALSXP>(weight_vec[index])){

            stop("RJMCMC_SMC.cpp. weight_vec[index] is NaN. probably from before");

          }

          if(Rcpp::traits::is_infinite<REALSXP>(weight_vec[index])){

            stop("RJMCMC_SMC.cpp. weight_vec[index] is infinite. probably from before");

          }
          // END CODE TO REMOVE

          // START CODE TO REMOVE
          if(NumericVector::is_na(out_weight)){

            stop("RJMCMC_SMC.cpp. out_weight is na");

          }

          if(Rcpp::traits::is_nan<REALSXP>(out_weight)){

            stop("RJMCMC_SMC.cpp. out_weight is NaN");

          }

          if(Rcpp::traits::is_infinite<REALSXP>(out_weight)){

            stop("RJMCMC_SMC.cpp. out_weight is infinite");

          }
          // END CODE TO REMOVE

          double old_weight = weight_vec[index] ;
          weight_vec[index] = weight_vec[index] * std::exp(out_weight) ;


          // START CODE TO REMOVE
          if(NumericVector::is_na(weight_vec[index])){

            Rcout << "this is weight_vec[index]   " << weight_vec[index] ;
            stop("RJMCMC_SMC.cpp. weight_vec[index] is na");

          }

          if(Rcpp::traits::is_nan<REALSXP>(weight_vec[index])){

            Rcout << "this is weight_vec[index]   " << weight_vec[index] ;

            stop("RJMCMC_SMC.cpp. weight_vec[index] is NaN");

          }

          if(Rcpp::traits::is_infinite<REALSXP>(weight_vec[index])){

            Rcout << "this is weight_vec[index]   " << weight_vec[index] << "\n" ;
            Rcout << "this is old_weight   " << old_weight << "\n" ;
            Rcout << "this is out_weight   " << out_weight << "\n" ;
            Rcout << "this is std::exp(out_weight)   " << std::exp(out_weight) << "\n" ;

            stop("RJMCMC_SMC.cpp. weight_vec[index] is infinite");

          }
          // END CODE TO REMOVE


          //store the new particles

          Svec[index] = S ;


          // Bvec has S+1 elements (S segments need S+1 breakpoints)
          // Bvec_final should also have S+1 elements
          NumericVector Bvec_final(S+1) ;
          IntegerVector Vvec_final(S) ;
          IntegerVector Zvec_final(S) ;
          IntegerVector Qvec_final(S) ;
          IntegerVector Fvec_final(S) ;

          // Bvec_final has the last breakpoint observed inside the previous Update Interval + all
          // sampled breakpoints inside the current  Update Interval

          double current_B_last = B_last[index] ;
          
          Bvec_final[0] = current_B_last ;
          Vvec_final[0] = Vvec[0] ;
          Zvec_final[0] = Zvec[0] ;
          Qvec_final[0] = Qvec[0] ;
          Fvec_final[0] = Fvec[0] ;
          
          if(S >= 2){

            for(int j=1; j<S;j++){

              Bvec_final[j] = Bvec[j] ;
              Vvec_final[j] = Vvec[j] ;
              Zvec_final[j] = Zvec[j] ;
              Qvec_final[j] = Qvec[j] ;
              Fvec_final[j] = Fvec[j] ;

            }

            // Copy the last breakpoint Bvec[S] to Bvec_final[S]
            // This is critical: Bvec[S] is the last breakpoint (usually > end_point for open segments)
            Bvec_final[S] = Bvec[S] ;

            // Update B_last: the last changepoint observed WITHIN the Update Interval

            double new_B_last = Bvec[S-1] ;

            if (new_B_last > end_point){
              Rcpp::stop("RJMCMC_SMC.cpp at end of loop: new_B_last = %f must be <= end_point = %f for particle %d", 
                         new_B_last, end_point, index);
            }
            B_last[index] = new_B_last ;

            //update the V state of the last complete segment (start and end of the seg inside the Update Interval)
            V_last_complete[index] = Vvec[S - 2] ;

          }else{

            // Case when S == 1: Bvec has 2 elements (Bvec[0] and Bvec[1])
            // We've already set Bvec_final[0], now set Bvec_final[1]
            if(S == 1){

              Bvec_final[1] = Bvec[1] ;

              // Update B_last: when S==1, Bvec[0] is t_star and Bvec[1] is the last breakpoint
              //B_last[index] = Bvec[1] ;

              B_last[index] = end_point-0.001;
            
            }
          
          }

          if (Bvec_final[0] > start_point){

            Rcout << "this is current_B_last   " << current_B_last << "\n" ;
            Rcout << "this is start_point   " << start_point << "\n" ;
            Rcout << "this is Bvec_final   " << Bvec_final << "\n" ;
            Rcout << "this is S   " << S << "\n" ;
            Rcout << "this is Bvec   " << Bvec << "\n" ;
            
            Rcpp::stop("RJMCMC_SMC.cpp at end of loop: Bvec_final[0] = %f must be <= start_point = %f for particle %d", 
                       Bvec_final[0], start_point, index);
          }


          //this copy is by value
          container_B[index] = clone(Bvec_final) ;
          container_V[index] = clone(Vvec_final) ;
          container_Z[index] = clone(Zvec_final) ;
          container_Q[index] = clone(Qvec_final) ;
          container_F[index] = clone(Fvec_final) ;
          


          // update the particle_index iterator
          particle_index = particle_index + 1 ;
          particle_count = particle_count - 1 ;


        }else{

          // DO NOTHING. No particles to update
        }


        if(particle_count == 0){

          break ;

        }

      }

    } //end of the  RJMCMC iteration

    if (particle_count > 0){
      Rcout << "MCMC did not update all the particles, please increase n_ite parameter " << "\n" ;
      Rcpp::stop("RJMCMC_SMC.cpp: MCMC did not update all the particles, please increase n_ite parameter ");
    }

    for(int i=0; i<particle_index_vec_local.length(); i++){

      int ii = particle_index_vec_local[i];

      if(B_last[ii] > end_point){

        Rcout << "this is B_last   " << B_last << "\n" ;

        Rcpp::stop("RJMCMC_SMC.cpp at 602: B_last[%d] = %f must be <= end_point = %f", 
                   ii, B_last[ii], end_point);
      }
    }

   // Validate breakpoint constraints for all Bvec in container_B
    for(int i = 0; i < particle_index_vec_local.length(); i++){
    
      int ii = particle_index_vec_local[i];
        
      NumericVector Bvec_check = container_B[ii];
      int S_check = Svec[ii];
      
      // Bvec should have S+1 elements
      if(Bvec_check.length() != S_check + 1){
        Rcpp::stop("RJMCMC_SMC.cpp: Bvec length mismatch. Expected %d, got %d for particle %d", 
                   S_check + 1, Bvec_check.length(), ii);
      }
      
      // Check 1: Bvec[0] must be < start_point (carry-over from previous UI)
      if( (S_check>1) & (Bvec_check[1] < start_point)){
        Rcpp::stop("RJMCMC_SMC.cpp: Bvec[1] = %f must be >= start_point = %f for particle %d", 
                   Bvec_check[1], start_point, ii);
      }

      if (Bvec_check[0] > start_point){
        Rcpp::stop("RJMCMC_SMC.cpp: Bvec[0] = %f must be <= start_point = %f for particle %d", 
                   Bvec_check[0], start_point, ii);
      }
      // Check 2: Last element (Bvec[S]) must be > end_point (open segment boundary)
      int last_idx = Bvec_check.length() - 1;  // This is S (since Bvec has S+1 elements)
      if(Bvec_check[last_idx] <= end_point){
        Rcout << "start_point = " << start_point << "\n";
        Rcout << "end_point = " << end_point << "\n";
        Rcout << "this is T_seg --> " << T_seg << "\n" ;
        Rcout << "this is Bvec --> " << Bvec_check << "\n" ;
        Rcpp::stop("RJMCMC_SMC.cpp: Last breakpoint Bvec[%d] = %f must be > end_point = %f for particle %d", 
                   last_idx, Bvec_check[last_idx], end_point, ii);
      }
      
      // Check 3: Elements from index 1 to S-1 must be between start_point and end_point (inclusive)
      for(int j = 1; j < S_check; j++){
        if(Bvec_check[j] < start_point || Bvec_check[j] > end_point){
          Rcpp::stop("RJMCMC_SMC.cpp: Bvec[%d] = %f must be between start_point = %f and end_point = %f for particle %d", 
                     j, Bvec_check[j], start_point, end_point, ii);
        }
      }
    }

    if(particle_count >0){

      Rcout<< "MCMC did not update all the particles, please increase n_ite parameter " << "\n" ;

    }

    // calling the R function breakpoints_sampling
    Function breakpoints_sampling("breakpoints_sampling");

    // select Breakpoints vector only for the particles updated in this routine

    List selected_breakpoint_vector = container_B[particle_index_vec] ;

    // START CODE TO DELETE
    for ( int j=0; j<selected_breakpoint_vector.length();j++){

      NumericVector x = selected_breakpoint_vector[j];

      for (int i = 0; i < x.size(); i++) {
        if (NumericVector::is_na(x[i])) {
          stop("RJSMCMC_SMC.cpp --> NaN value detected in the vector.");
        }
      }
    }
    //END CODE TO DELETE

    int sample_size = 5 ;

    List generated_breakpoint_list = breakpoints_sampling(
                         Named("start_point") = start_point,
                         Named("end_point") = end_point,
                         Named("breakpoint_list") = selected_breakpoint_vector,
                         Named("sample_size") = sample_size) ;

    // compute correction weight
    double weigths_corr_value = weigths_correction(
                       T_seg,
                       generated_breakpoint_list,
                       sample_size,
                       K,
                       keyvec,
                       etavec,
                       key0vec,
                       eta0vec,
                       probvec_Z,
                       probvec_F,
                       P0) ;

    // correct weights for the particles updated in this routine

    double index_value ;

    for(int i=0 ; i<particle_index_vec.length(); i++){

      index_value = particle_index_vec[i] ;

      weight_vec[index_value] = weight_vec[index_value] * std::exp(weigths_corr_value);

      if(NumericVector::is_na(weight_vec[index_value])){

        
        Rcout << "this is index_value   " << index_value << "\n"  ;
        Rcout << "this is weight_vec[index_value]   " << weight_vec[index_value] << "\n"  ;

        stop("RJMCMC_SMC.cpp. weight is na line 721");

        }

      if(NumericVector::is_na(weigths_corr_value)){

        stop("RJMCMC_SMC.cpp. weigths_corr_value is na");

      }

      if(Rcpp::traits::is_nan<REALSXP>(weigths_corr_value)){

        stop("RJMCMC_SMC.cpp. weigths_corr_value is NaN");

      }

      if(Rcpp::traits::is_infinite<REALSXP>(weigths_corr_value)){

        stop("RJMCMC_SMC.cpp. weigths_corr_value is infinite");

      }

    }

}









