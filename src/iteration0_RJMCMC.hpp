#ifndef iteration0_header
#define iteration0_header

using namespace Rcpp;

List iteration0_RJMCMC (   const NumericVector& T_seg,
                           const IntegerVector& Y_seg,
                           const int& minimum_n,
                           const double& start_point,
                           const double& end_point,
                           const double& t_star,
                           const int& K,
                           const int& W,
                           const int& U,
                           const bool& empty_mix,
                           const NumericVector& probvec_V,
                           const NumericVector& probvec_Z,
                           const NumericVector& probvec_Q,
                           const NumericVector& probvec_F,
                           const NumericVector& alphavec,
                           const NumericVector& muvec,
                           const NumericVector& keyvec,
                           const NumericVector& etavec,
                           const NumericVector& key0vec,
                           const NumericVector& eta0vec,
                           const NumericMatrix lambdamat,
                           const double& P0,
                           const int& num_logs,
                           const double& max_range,
                           const int& Smax);


#endif
