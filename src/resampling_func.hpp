#ifndef resampling_func_header
#define resampling_func_header


using namespace Rcpp;


List resampling_func(const NumericVector& weight_vec,
                     const List& container_B,
                     const List& container_V,
                     const List& container_Z,
                     const List& container_Q,
                     const List& container_F,
                     const IntegerVector& Svec,
                     const IntegerVector& V_last_complete,
                     const NumericVector& B_last,
                     int n_particle);

#endif
