#include <Rcpp.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <algorithm>
#include "toolbox1.hpp"
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//' @title Function for approximating the log probability of observing "n" obs in semi-segment [tau0,split_U) given
//' the total length [tau0,tau1)
//' @param count (double& int) number of observations within the selected segment
//' @param split_U (const& double) point where the update interval ends
//' @param alpha (const& double) shape parameter for the rate
//' @param mu (const& double) mean parameter for the rate
//' @param minimum_n (const& int) minimum number of observations allowed in a non-empty segment
//' @param tau0 (const& double) left bound of the segment
//' @param tau1 (const& double) right bound of the segment
//' @param error (const& double) approximation error
//' @return the value of the log density
// [[Rcpp::export]]
double prob_approx( const double& count,
                    const double& split_U,
                    const double& alpha,
                    const double& mu,
                    const int& minimum_n,
                    const double& tau0,
                    const double& tau1,
                    const double& error
){



  //initialize objects

  double log_first_side;

  int max_value = std::max(minimum_n,(int)count)   ;

  double old_value = -1 ;
  double new_value = 0  ;

  double  log_second_side=0,second_side=0  ;

  //initialize step_vector for approximating the infinite sum

  int step = 1000 ; int step_index = 1 ;

  IntegerVector m_step = seq(0,step) ;

  log_first_side = count * ( std::log(split_U- tau0) - std::log(tau1-split_U) )

                 + alpha * (std::log(alpha/mu)  - std::log(alpha/mu + tau1 - tau0)  )

                 + minimum_n * (std::log(alpha/mu + tau1 - tau0) -std::log(tau1 - tau0) ) ;



  // stop the sum when the difference between the sum at previous iteration snd the sum at current iteration is less than "error"
  while( abs(new_value-old_value) >error){

      log_second_side=0;
      second_side=0  ;

    old_value = new_value ;

    //in order to compute the second side I need to compute the binomial coefficient
    // binom(a+max_value-delta+m_step[i]-1,max_value-delta+m_step[i]), which changes at each iteration

    double log_binom_coeff = std::log(alpha) ;

    for(int i=1; i<max_value-minimum_n;i++){

      log_binom_coeff += std::log(alpha+i) ;

    }

    int updating_value = max_value - minimum_n ;

    //compute the second part of the density: it is composed by two binomial coefficient, 1 involving only interger and
    //one involving alpha which is a real number (hence "binomial_coefficient" function from boost cannot be used)

    for(int i=0; i<m_step.length();i++){



      log_second_side  =  log_binom_coeff -  lgamma(updating_value+1)

                       + lgamma(m_step[i] + max_value+1) - lgamma(count+1) - lgamma(m_step[i] + max_value - count +1)

                       + ((double)m_step[i] + (double)max_value) *std::log( (tau1 - split_U)/(alpha/mu + tau1 - tau0))

                       + log_first_side ;

     if(log_second_side<-700.0){log_second_side = -700.0 ; }

      second_side+= std::exp(log_second_side) ;



    //update log_binom_coeff

    updating_value = max_value - minimum_n + (m_step[i]+1) ;

    log_binom_coeff += std::log(alpha + updating_value -1) ;

    }

    if(second_side==R_NegInf){second_side =0.0 ;}



    new_value =   second_side  ;


    //new_value =  first_side * second_side  ;

    if(new_value>1){

      new_value=0.99 ;

    }

    step_index = step_index + 1 ;

    m_step = seq(0,step * step_index) ;


  }

  if(R_isnancpp(new_value)){

    Rcout << "log_first_side " << log_first_side << "\n" ;
    Rcout << "second_side " << second_side << "\n" ;
    Rcout << "count " << count << "\n" ;
    Rcout << "split_u " << split_U << "\n" ;
    Rcout << "alpha " << alpha << "\n" ;
    Rcout << "mu " << mu << "\n" ;
    Rcout << "tau0 " << tau0 << "\n" ;
    Rcout << "tau1 " << tau1 << "\n" ;


    stop("returned NaN when approximating the density for the count");}

  if(NumericVector::is_na(new_value)){stop("returned NA/NaN");}

  if(Rcpp::traits::is_infinite<REALSXP>(new_value)){

    Rcout << "log_first_side " << log_first_side << "\n" ;
    Rcout << "second_side " << second_side << "\n" ;
    Rcout << "count " << count << "\n" ;
    Rcout << "split_u " << split_U << "\n" ;
    Rcout << "alpha " << alpha << "\n" ;
    Rcout << "mu " << mu << "\n" ;
    Rcout << "tau0 " << tau0 << "\n" ;
    Rcout << "tau1 " << tau1 << "\n" ;

    stop("returned inf");}


  return std::log(new_value) ;

}
