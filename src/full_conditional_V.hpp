#ifndef full_conditional_V_header
#define full_conditional_V_header

using namespace Rcpp;


List full_conditional_V(const int& V,
                        const int& U,
                        const int& num_logs,
                        const NumericVector& T_seg,
                        const IntegerVector& Y_seg,
                        const double& LB,
                        const double& UB,
                        const NumericMatrix& lambdamat,
                        const NumericVector& probvec_V,
                        int V_left,
                        const double& P0,
                        const double& end_point=0.0,
                        const bool open_segment=false,
                        const bool sample_V=false);

#endif

