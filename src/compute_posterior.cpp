#include <Rcpp.h>
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

  NumericMatrix posterior_matrix(num_discr_intervals, num_states);

  for(int i=0; i<num_particles ; i++){

    for(int j=0; j<num_discr_intervals ; j++){

      posterior_matrix(j,state_container(j,i)) += weight_vec[i] ;

      }

  }


  return posterior_matrix;

}

