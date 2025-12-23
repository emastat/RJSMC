#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
//'
//' Function for updating the posterior distribution of state variables
//' @param num_discr_intervals (int) number of discretization intervals
//' @param num_particles (int) number of particles
//' @param state_container (IntegerMatrix) with inferred state values for each discretization interval, for each particle
//' @param num_states (int) number of states
//' @param weight_vec (NumericVector) weights assigned to each particle
//' @export
// [[Rcpp::export]]

NumericMatrix compute_posterior(const int num_discr_intervals,
                                const int num_particles,
                                const IntegerMatrix state_container,
                                const int num_states,
                                const NumericVector weight_vec) {

  // VALIDATION CHECK 1: Verify weights sum to 1 (with tolerance for floating point)
  // Threshold for accepting very small negative values (should be rounded to 0 before reaching here)
  const double weight_threshold = 1e-5;
  double sum_weights = 0.0;
  for(int i = 0; i < num_particles; i++){
    if(NumericVector::is_na(weight_vec[i]) || Rcpp::traits::is_nan<REALSXP>(weight_vec[i])){
      Rcpp::stop("compute_posterior: weight_vec contains NA or NaN at index %d", i);
    }
    // Only error if negative value is significantly negative (beyond numerical precision)
    if(weight_vec[i] < -weight_threshold){
      Rcpp::stop("compute_posterior: weight_vec contains negative value at index %d: %.15f (threshold: %.1e)", 
                 i, weight_vec[i], weight_threshold);
    }
    // Treat very small negative values as 0 when summing (shouldn't happen if normalization is correct)
    double weight_value = (weight_vec[i] < 0.0) ? 0.0 : weight_vec[i];
    sum_weights += weight_value;
  }
  
  // Tolerance for sum check: allow small numerical approximation errors
  const double tolerance = 1e-5;
  if(std::abs(sum_weights - 1.0) > tolerance){
    Rcpp::stop("compute_posterior: weights do not sum to 1. Sum=%.15f, expected=1.0, difference=%.15f (tolerance: %.1e)", 
               sum_weights, std::abs(sum_weights - 1.0), tolerance);
  }

  // VALIDATION CHECK 2: Verify state_container dimensions match expected
  if(state_container.nrow() != num_discr_intervals){
    Rcpp::stop("compute_posterior: state_container.nrow()=%d does not match num_discr_intervals=%d", 
               state_container.nrow(), num_discr_intervals);
  }
  if(state_container.ncol() != num_particles){
    Rcpp::stop("compute_posterior: state_container.ncol()=%d does not match num_particles=%d", 
               state_container.ncol(), num_particles);
  }

  NumericMatrix posterior_matrix(num_discr_intervals, num_states);

  for(int i=0; i<num_particles ; i++){

    for(int j=0; j<num_discr_intervals ; j++){

      // VALIDATION CHECK 3: Verify state index is valid
      int state_idx = state_container(j,i);
      if(state_idx < 0 || state_idx >= num_states){
        Rcpp::stop("compute_posterior: Invalid state index %d at position [%d,%d]. Valid range is [0, %d)", 
                   state_idx, j, i, num_states);
      }

      posterior_matrix(j,state_idx) += weight_vec[i] ;

      }

  }

  // VALIDATION CHECK 4: Verify each row sums correctly (should equal 1.0 within tolerance)
  // Use same tolerance as sum check (1e-5) to allow for numerical approximation errors
  const double row_tolerance = 1e-5;
  for(int j = 0; j < num_discr_intervals; j++){
    double row_sum = 0.0;
    for(int s = 0; s < num_states; s++){
      row_sum += posterior_matrix(j, s);
    }
    if(std::abs(row_sum - 1.0) > row_tolerance){
      Rcpp::stop("compute_posterior: Row %d does not sum to 1.0. Sum=%.15f, difference=%.15f (tolerance: %.1e)", 
                 j, row_sum, std::abs(row_sum - 1.0), row_tolerance);
    }
  }

  return posterior_matrix;

}

