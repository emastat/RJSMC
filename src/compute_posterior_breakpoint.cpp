#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//'
//' Function for computing the posterior distribution of breakpoint locations
//' @param num_discr_intervals (int) number of discretization intervals
//' @param num_particles (int) number of particles
//' @param state_container_B (NumericMatrix) matrix indicating breakpoint presence: 1 if particle j has a breakpoint in interval i, 0 otherwise
//' @param weight_vec (NumericVector) weights assigned to each particle (should sum to 1)
//' @return NumericVector posterior probability that a breakpoint occurs in each discretization interval
//' @export
// [[Rcpp::export]]

NumericVector compute_posterior_breakpoint(const int num_discr_intervals,
                                          const int num_particles,
                                          const NumericMatrix state_container_B,
                                          const NumericVector weight_vec) {

  // VALIDATION CHECK 1: Verify weights sum to 1 (with tolerance for floating point)
  // Threshold for accepting very small negative values (should be rounded to 0 before reaching here)
  const double weight_threshold = 1e-5;
  double sum_weights = 0.0;
  for(int i = 0; i < num_particles; i++){
    if(NumericVector::is_na(weight_vec[i]) || Rcpp::traits::is_nan<REALSXP>(weight_vec[i])){
      Rcpp::stop("compute_posterior_breakpoint: weight_vec contains NA or NaN at index %d", i);
    }
    // Only error if negative value is significantly negative (beyond numerical precision)
    if(weight_vec[i] < -weight_threshold){
      Rcpp::stop("compute_posterior_breakpoint: weight_vec contains negative value at index %d: %.15f (threshold: %.1e)", 
                 i, weight_vec[i], weight_threshold);
    }
    // Treat very small negative values as 0 when summing (shouldn't happen if normalization is correct)
    double weight_value = (weight_vec[i] < 0.0) ? 0.0 : weight_vec[i];
    sum_weights += weight_value;
  }
  
  // Tolerance for sum check: allow small numerical approximation errors
  const double tolerance = 1e-5;
  if(std::abs(sum_weights - 1.0) > tolerance){
    Rcpp::stop("compute_posterior_breakpoint: weights do not sum to 1. Sum=%.15f, expected=1.0, difference=%.15f (tolerance: %.1e)", 
               sum_weights, std::abs(sum_weights - 1.0), tolerance);
  }

  // VALIDATION CHECK 2: Verify state_container_B dimensions match expected
  if(state_container_B.nrow() != num_discr_intervals){
    Rcpp::stop("compute_posterior_breakpoint: state_container_B.nrow()=%d does not match num_discr_intervals=%d", 
               state_container_B.nrow(), num_discr_intervals);
  }
  if(state_container_B.ncol() != num_particles){
    Rcpp::stop("compute_posterior_breakpoint: state_container_B.ncol()=%d does not match num_particles=%d", 
               state_container_B.ncol(), num_particles);
  }

  // Initialize posterior vector with zeros
  NumericVector posterior_B_vec(num_discr_intervals);
  for(int i = 0; i < num_discr_intervals; i++){
    posterior_B_vec[i] = 0.0;
  }

  // Compute posterior distribution: for each interval, sum weights of particles with breakpoints
  for(int i = 0; i < num_discr_intervals; i++){
    for(int j = 0; j < num_particles; j++){
      // Check if particle j has a breakpoint in interval i
      if(state_container_B(i, j) == 1.0){
        // Add weight of particle j to the posterior for interval i
        double weight_value = (weight_vec[j] < 0.0) ? 0.0 : weight_vec[j];
        posterior_B_vec[i] += weight_value;
      }
    }
  }

  return posterior_B_vec;

}

