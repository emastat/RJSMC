#include <algorithm>
#include <Rcpp.h>
#include "extract_index.hpp"
#include "open_negative_binomial.hpp"
#include "toolbox2.hpp"

using namespace Rcpp;



// [[Rcpp::plugins(cpp11)]]

//' @title  posterior_evaluation
//' Function for evaluating the joint posterior distribution of all the quantities related to a closed or open  segment (T,Y,V,Z,Q,F)
//'
//' @param T_seg Vector with time stamps within the update interval
//' @param Y_seg Vector with the messages within the update interval
//' @param V (const int&) Value of the current V if present; if not set it 0.
//' @param Z (const int&) Value of the current Z if present; if not set it 0.
//' @param Q (const int&) Value of the current Q if present; if not set it 0.
//' @param F (const int&) Value of the current F if present; if not set it 0.
//' @param U (const int&) Number of levels for the V state
//' @param K (const int&) Number of levels of the Z state
//' @param W (const int&) Number of levels for the Q state
//' @param LB (const double&) lower bound of the selected segment
//' @param UB (const double&) upper bound of the selected segment
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
//' @param V_left (const int&) V value for the previous segment
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param minimum_n (const int&) minimum number of observations that must be observed in a non-empty segment
//' @param count (const double&) number of observarions inside the selected segment
//' @param open_segment (const bool, default is false) if false  the segment is considered totally observed. If true the semgent is considered partially observed; this affects the likelihood computation
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the stops to be observed
//' @return  List with 6 elements
//' \describe{
//' \item{ldens_tmsp}{double,likelihood of the time-stamps positions within the selected segment }
//' \item{ldens_count}{double, likelihood for the count of observations withing the selected segment }
//' \item{ldens_event}{double,  likelihood for the messages within the selected segment}
//' \item{ldens_empty}{double, prior probability to observe state F in the selected segment}
//' \item{ldens_L}{double, prior for  the length of the selected segment}
//' \item{ldens_dynamics}{double, prior probability to observe state V, Z and Q in the selected segment}
//' }
//' @export
// [[Rcpp::export]]
double posterior_evaluation(const NumericVector& T_seg,
                            const IntegerVector& Y_seg,
                            const int& V,
                            const int& Z,
                            const int& Q,
                            const int& F,
                            const int& U,
                            const int& K,
                            const int& W,
                            const double& LB,
                            const double& UB,
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
                            const int& V_left,
                            const double& P0,
                            const int& minimum_n,
                            const double& count,
                            const double& end_point=0.0,
                            const bool open_segment=false){



  // L: length of the selected segment
  double L = UB - LB ;

  //TIME-STAMPS POSITION: T_s1... Ts,n_s
  double ldens_tmsp =std::lgamma(count + 1) - count * std::log(L)  ;

  //STATE_DYNAMICS:  V, Z, Q, F

  double ldens_dynamics;


  //initialize the density values
  double ldens_count  ,  ldens_empty, ldens_event , ldens_L ;

  if(count>0){


    if(open_segment==false){

    // ldens_count: (collapsed) log density for the count of obs inside a CLOSED segment
    ldens_count = R::dnbinom(count-minimum_n,
                             alphavec[Q-1],
                             (alphavec[Q-1]/muvec[Q-1])/(alphavec[Q-1]/muvec[Q-1]+ L)
                             ,true) ;

    }else{

      // ldens_count: (collapsed) log density for the count of obs inside an OPEN segment
      ldens_count = prob_approx(count,
                                end_point,
                                alphavec[Q-1],
                                muvec[Q-1],
                                minimum_n,
                                LB,
                                UB,
                                0.0001);

    }

    //initialize out_extract: vector with indexes of the observations contained in the segment
    NumericVector out_extract(3) ;

    // compute log_vector: vector with messages whithin the selected segment

    if(open_segment==false){

      out_extract =clone( extract_index(T_seg,LB,UB)) ; // observations found in the whole segment

    }else{

       out_extract =clone( extract_index(T_seg,LB,end_point)) ; // observations found in the partial segment
    }

    IntegerVector log_vector = Y_seg[Range(out_extract[1],out_extract[2])] ;

    // frequency table for the messages
    IntegerVector dj_vec = table_cpp(log_vector,num_logs,false) ;


    //ldens_event: likelihood for the messages included in the selected segment
    for(int j=0 ; j < num_logs ; j++){

      ldens_event = ldens_event + std::log(lambdamat(V-1,j) ) * (double)dj_vec[j] ;

    }

    //ldens_empty: prior density for observing an empty segment (is 0 because we are in non-empty segment)
    ldens_empty = 0 ;


    //ldens_dynamics: prior density for the state V,Z,Q

    if(V_left==0){

      ldens_dynamics = std::log(probvec_V[V-1]) + std::log(probvec_Z[Z-1]) + std::log(probvec_Q[Q-1]) ;

      if(Rcpp::traits::is_infinite<REALSXP>(ldens_dynamics)){

        double prob_V = std::log(probvec_V[V-1]) ;
        double prob_Z = std::log(probvec_Z[Z-1]) ;
        double prob_Q = std::log(probvec_Q[Q-1]) ;

        Rcout << "prob_V: " << prob_V << "\n" ;
        Rcout << "prob_Z: " << prob_Z << "\n" ;
        Rcout << "prob_Q: " << prob_Q << "\n" ;
        Rcout << "V: " << V << "\n" ;
        Rcout << "Z: " << Z << "\n" ;
        Rcout << "Q: " << Q << "\n" ;

        stop("ldens_dynamic is infinite: posterior_evluaztion.cpp");

      }

    }else{

      ldens_dynamics = std::log(1-P0) + std::log(probvec_V[V-1]) + std::log(probvec_Z[Z-1]) + std::log(probvec_Q[Q-1]) ;

      if(Rcpp::traits::is_infinite<REALSXP>(ldens_dynamics)){

        double prob_V = std::log(probvec_V[V-1]) ;
        double prob_Z = std::log(probvec_Z[Z-1]) ;
        double prob_Q = std::log(probvec_Q[Q-1]) ;
        double prob_P0 = std::log(1-P0) ;

        Rcout << "prob_V: " << prob_V << "\n" ;
        Rcout << "prob_Z: " << prob_Z << "\n" ;
        Rcout << "prob_Q: " << prob_Q << "\n" ;
        Rcout << "prob_P0: " << prob_P0 << "\n" ;
        Rcout << "V: " << V << "\n" ;
        Rcout << "Z: " << Z << "\n" ;
        Rcout << "Q: " << Q << "\n" ;

        stop("ldens_dynamic is infinite: posterior_evluaztion.cpp");
      }

    }

  }else{

      ldens_count = 0  ;

      ldens_event = 0 ;

      ldens_empty = std::log(probvec_F[F-1]) ;

      ldens_dynamics = std::log(P0) ;

  }

    // ldens_L: prior density for the length of the selected segment

    if(count>0){

      //non-empty segment
      ldens_L = R::dgamma( L, keyvec[Z-1],etavec[ Z-1 ]/keyvec[Z-1],true );

    }else{

      //empty segment
      ldens_L = R::dgamma( L, key0vec[F-1] ,eta0vec[ F-1 ]/key0vec[F-1],true );   }




    if(NumericVector::is_na(ldens_tmsp)){
      Rcout << "ldens_tmsp: " << ldens_tmsp << "\n" ;
      stop("ldens_tmsp is na: posterior_evaluztion.cpp");}

    if(Rcpp::traits::is_infinite<REALSXP>(ldens_tmsp)){
      Rcout << "ldens_tmsp: " << ldens_tmsp << "\n" ;
      stop("ldens_tmsp is infinite :posterior_evaluztion.cpp");}


    if(NumericVector::is_na(ldens_count)){
      Rcout << "ldens_count: " << ldens_count << "\n" ;
      stop("ldens_count is na: posterior_evaluztion.cpp");}

    if(Rcpp::traits::is_infinite<REALSXP>(ldens_count)){
      Rcout << "ldens_count: " << ldens_count << "\n" ;
      stop("ldens_count is infinite: posterior_evaluztion.cpp");}



    if(NumericVector::is_na(ldens_event)){
      Rcout << "ldens_event: " << ldens_event << "\n" ;
      stop("ldens_event is na: posterior_evaluztion.cpp");}


    if(Rcpp::traits::is_infinite<REALSXP>(ldens_event)){
      Rcout << "ldens_event: " << ldens_event << "\n" ;
      stop("ldens_event is infinite: posterior_evaluztion.cpp");}


    if(NumericVector::is_na(ldens_empty)){
      Rcout << "ldens_empty: " << ldens_empty << "\n" ;
      stop("ldens_empty is na: posterior_evaluztion.cpp");}

    if(Rcpp::traits::is_infinite<REALSXP>(ldens_empty)){
      Rcout << "ldens_empty: " << ldens_empty << "\n" ;
      stop("ldens_empty is infinite: posterior_evaluztion.cpp");}


    if(NumericVector::is_na(ldens_L)){
      Rcout << "ldens_L: " << ldens_L << "\n" ;
      stop("ldens_L is na: posterior_evaluztion.cpp");}

    if(Rcpp::traits::is_infinite<REALSXP>(ldens_L)){
      Rcout << "ldens_L: " << ldens_L << "\n" ;
      stop("ldens_L is infinite: posterior_evaluztion.cpp");}


    if(NumericVector::is_na(ldens_dynamics)){
      Rcout << "ldens_dynamics: " << ldens_dynamics << "\n" ;
      stop("ldens_dynamics is na: posterior_evaluztion.cpp");}

    if(Rcpp::traits::is_infinite<REALSXP>(ldens_dynamics)){
      Rcout << "ldens_dynamics: " << ldens_dynamics << "\n" ;
      stop("ldens_dynamics is infinite: posterior_evaluztion.cpp");}


  return(ldens_tmsp +
         ldens_count +
         ldens_event +
         ldens_empty +
         ldens_L +
         ldens_dynamics ) ;



}











