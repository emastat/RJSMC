#ifndef collapsed_full_conditional_Q_header
#define collapsed_full_conditional_Q_header

using namespace Rcpp;



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
                        const bool sample_state=false);


#endif