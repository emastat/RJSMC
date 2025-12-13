#ifndef incremental_weight_header
#define incremental_weight_header

using namespace Rcpp;


double incremental_weight(const NumericVector& T_seg,
                          const IntegerVector& Y_seg,
                          const int& V,
                          const int& Z,
                          const int& Q,
                          const int& F,
                          const int& U,
                          const int& K,
                          const int& W,
                          const double& B_last,
                          const double& B1,
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
                          const int& V_last_complete,
                          const double& P0,
                          const int& minimum_n,
                          const double& t_star,
                          const double& end_point=0.0,
                          const bool open_segment=false);


#endif
  
  
