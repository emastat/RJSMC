#ifndef initial_V_header
#define initial_V_header

using namespace Rcpp;

List initial_state ( const NumericVector& T_seg,
                     const IntegerVector& Y_seg,  
                     const int& S,
                     const NumericVector& Bvec,
                     const double& end_point,  
                     const NumericVector& probvec_V,
                     const NumericVector& probvec_Q,
                     const NumericVector& probvec_Z,
                     const NumericVector& probvec_F,
                     const NumericMatrix& lambdamat,
                     const int& U,
                     const int& W,
                     const int& K,
                     const NumericVector& alphavec,
                     const NumericVector& muvec,
                     const NumericVector& keyvec,
                     const NumericVector& etavec,
                     const NumericVector& key0vec,
                     const NumericVector& eta0vec,
                     const double& P0,
                     const int& num_logs,
                     const bool& empty_mix,
                     const int& minimum_n);


#endif

