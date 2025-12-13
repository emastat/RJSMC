#include <Rcpp.h>
#include "toolbox2.hpp"
#include "open_negative_binomial.hpp"
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
//' @title collapsed_full_conditional_Q
//' Function for sampling the Q  state for a  new non-empty segment from its  collpsed full conditional and/or evaluate it at the sampled point
//' It is "collapsed" because the likelihood part is not Gamma(M,alpha_q,alpha_q/mu_q) but
//' shifted Neg-Bin (N;alpha_q,L/(alpha_q/mu_q + L)). M is integrated out because not needed
//' @param Q (const int&) Value of the current Q if present; if not set it 0.
//' @param W (const int&) Number of levels for the Q state
//' @param count (const double&) number of observations within the selected segment
//' @param LB (const double&) lower bound of the selected segment
//' @param UB (const double&) upper bound of the selected segment
//' @param alphavec (const NumericVector&) Shape parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param muvec (const NumericVector&) Mean parameters ruling the rate of occurrence of segment for the different levels of the Q state
//' @param probvec_Q (const NumericVector&) Probability mass for the Q state distribution
//' @param minimum_n minimum number of observations that must be observed in a non-empty segment
//' @param open_segment (const bool) if false (by default) the segment is considered totally observed. If true the semgent is considered partially observed; this affects the likelihood computation
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the stops to be observed
//' @param sample_Q ( const bool=false) if true a new Q is sampled and the full conditional is evaluated at the new value. If false (by default) the full condititional is evaluated at the input Q
//' @return A list containting 2 elements:
//' \describe{
//' \item{Q_new}{int, the sample Q state for the selected segment}
//' \item{eval_densQ}{double, value of the log full conditional for either the new sampled Q or the current Q inputed ( depending on sample_Q)}
//' }
//' @export
// [[Rcpp::export]]
List collapsed_full_conditional_Q(const int& Q,
                                  const int& W,
                                  const double& count,
                                  const double& LB,
                                  const double& UB,
                                  const NumericVector& alphavec,
                                  const NumericVector& muvec,
                                  const NumericVector& probvec_Q,
                                  const int& minimum_n,
                                  const bool open_segment=false,
                                        double end_point= 0.0,
                                  const bool sample_Q=false){


  //L: length of the selected segment
  double L = UB - LB ;

  // initialize Q and other  required quantities

  int Q_new=0 ;

  double eval_densQ, log_neg_binom ;

  NumericVector log_densQ(W), densQ(W);

  // compute probability weights for the full conditional

  for(int i=0; i<W; i++){

    if(open_segment==false){

      //the segment if fully observed: the likelihood for the count is a shifted Neg-binomial
      log_neg_binom = R::dnbinom(count-minimum_n, alphavec[i], (alphavec[i]/muvec[i])/(alphavec[i]/muvec[i]+ L),true) ;

    }else{

      //the segment is observed only partially: the likelihood for the count must be approximated

      log_neg_binom = prob_approx(count,
                                  end_point,
                                  alphavec[i],
                                  muvec[i],
                                  minimum_n,
                                  LB,
                                  UB,
                                  0.0001);

    }


    log_densQ[i] = std::log( probvec_Q[i] ) + log_neg_binom  ;

  }

  // employ the log-sum-exp trick

  densQ = stabilize(log_densQ)  ; // final probability vector


  // a new Q must be sampled - log full conditional evaluated in this point
  if(sample_Q==true){

    // sample Q

    Q_new = as<int>( Rcpp::sample( W, 1, false, densQ, true ) ) ;

    // logarithm of the evaluated density at the sampled point
    eval_densQ = std::log(densQ[ Q_new-1 ]);

    // return a list with the new sample Q and the value of the full conditional in this value

    return( List::create(Named("Q")=Q_new,
                         Named("eval_densQ")=eval_densQ)) ;

  }else{

    // Q must NOT be sampled - log full conditional is evaluated at the current Q

    // logarithm of the evaluated the density in the current point

    eval_densQ = std::log(densQ[ Q-1 ]);


    // reuturn a list with the current Q and the value of the full conditional in this value

    return( List::create(Named("Q")=Q,
                         Named("eval_densQ")=eval_densQ)) ;


  }





}







