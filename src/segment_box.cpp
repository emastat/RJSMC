#include <Rcpp.h>
#include "full_conditional_V.hpp"
#include "full_conditional_Z.hpp"
#include "collapsed_full_conditional_Q.hpp"
#include "full_conditional_F.hpp"
#include "posterior_evaluation.hpp"
#include "extract_index.hpp"

using namespace Rcpp;

//' @title center_segment_box
//' Function for computing all the necessary quantities needed for updating/evaluating a CLOSED segment or an OPEN segment with at least 1 obs inside.
//' This function will either sample new states V,Z,Q,F or use the current ones to evaluate the joint proposal and
//'  and the joint posterior distribution for these values
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
//' @param open_segment (const bool, default is false) if false  the segment is considered totally observed. If true the semgent is considered partially observed; this affects the likelihood computation
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the stops to be observed
//' @param sample_state (const bool=false) if true new V,Z,Q,F are sampled and the full conditional are evaluated at the new values. If false (by default) the full condititionals are evaluated at the input provided
//' @return List with 6 elements
//' \describe{
//' \item{V}{int,Value of the current or new V}
//' \item{Z}{int,Value of the current or new Z}
//' \item{Q}{int,Value of the current or new Q}
//' \item{F}{int,Value of the current or new F}
//' \item{sum_proposals}{double,sum of the log proposal for the current or new values of V, Z, Q and F }
//' \item{out_posterior_evaluation}{double, joint log posterior distribuion for the quantities related to the selected segment }
//' }
//' @export
// [[Rcpp::export]]
List        segment_box(const NumericVector& T_seg,
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
                        const bool open_segment=false,
                        const double& end_point=0.0,
                        const bool sample_state=false){



  // L: length of the selected segment
   double L = UB - LB ;

  // extract number of observations within the selected segment
  NumericVector out_extract(3) ;

  if(open_segment==false){

       out_extract =clone( extract_index(T_seg,LB,UB)) ;

  }else{

       out_extract =clone( extract_index(T_seg,LB,end_point)) ;
  }

  double count = out_extract[0] ; // number of obs in [LB;UB/end_point)


  if(count>0){

    if(open_segment==false & count<minimum_n){

      stop("a closed segment must have at least minimum_n observations");

    }


    // sample V and/or evaluate the full conditional and the posterior distribution in V (either new or old)

    List out_V =  full_conditional_V(V,
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
                                     sample_state);

    //sample V and or evaluate the full conditional and the posterior distribution in V (either new or old)

    List out_Z = full_conditional_Z(Z,
                                    K,
                                    L,
                                    keyvec,
                                    etavec,
                                    probvec_Z,
                                    sample_state) ;

    //sample Q and or evaluate the  collapsed full conditional and the posterior distribution in Q (either new or old)

    List out_Q = collapsed_full_conditional_Q(Q,
                                              W,
                                              count,
                                              LB,
                                              UB,
                                              alphavec,
                                              muvec,
                                              probvec_Q,
                                              minimum_n,
                                              open_segment,
                                              end_point,
                                              sample_state);



    // sum all the log proposal distributions

    double sum_proposals = as<double>(out_V["eval_densV"])

                         + as<double>(out_Z["eval_densZ"])

                         + as<double>(out_Q["eval_densQ"]) ;


    // evaluate the posterior distribution for the quantities related to the selected segment


    double  out_post_evaluation = posterior_evaluation(
                                                T_seg,
                                                Y_seg,
                                                as<int>(out_V["V"]),
                                                as<int>(out_Z["Z"]),
                                                as<int>(out_Q["Q"]),
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
                                                count,
                                                end_point,
                                                open_segment);

    if(NumericVector::is_na(out_post_evaluation)){stop("returned NA/NaN in poterior evaluation-  segment box row 197");}
    if(Rcpp::traits::is_infinite<REALSXP>(out_post_evaluation)){stop("returned inf in poterior evaluation-  segment box row 197");}



    return( List::create(Named("V")=out_V["V"],
                         Named("Z")=out_Z["Z"],
                         Named("Q")=out_Q["Q"],
                         Named("F")=0,
                         Named("sum_proposals")=sum_proposals,
                         Named("sum_posterior")=out_post_evaluation)) ;


  }else{

   //sample F and or evaluate the full conditional and the posterior distribution in F (either new or old)

    List out_F = full_conditional_F(F,
                                    L,
                                    key0vec,
                                    eta0vec,
                                    probvec_F,
                                    sample_state);


  // sum all the log proposal distributions (only eval_densF because the segment is empty)

  double sum_proposals = as<double>(out_F["eval_densF"]) ;

  // evaluate the posterior distribuion for the quantities related to the selected segment

  double out_post_evaluation = posterior_evaluation(T_seg,
                                              Y_seg,
                                              0,
                                              0,
                                              0,
                                              as<int>(out_F["F"]),
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
                                              count,
                                              end_point,
                                              open_segment);


  return( List::create(Named("V")=0,
                       Named("Z")=0,
                       Named("Q")=0,
                       Named("F")=out_F["F"],
                       Named("sum_proposals")=sum_proposals,
                       Named("sum_posterior")=out_post_evaluation)) ;

  }

}










