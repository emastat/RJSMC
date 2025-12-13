#ifndef new_forward_function_cpp_header
#define new_forward_function_cpp_header

using namespace Rcpp;

void forward_function(      int& S, 
                            const NumericVector& Tvec, 
                            const IntegerVector& Yvec, 
                            const int& U,
                            const int& K,
                            const int& W,
                            const NumericVector& probvec_V,
                            const NumericVector& probvec_Z,
                            const NumericVector& probvec_Q,
                            const NumericVector& probvec_F,
                            const NumericMatrix& lambdamat,
                            const NumericVector& keyvec,
                            const NumericVector& etavec,
                            const NumericVector& key0vec,
                            const NumericVector& eta0vec,
                            const NumericVector& alphavec,
                            const NumericVector& muvec,
                            const int& num_logs, 
                            NumericVector& Bvec, 
                            IntegerVector& Zvec, 
                            IntegerVector& Qvec,
                            IntegerVector& Vvec,
                            IntegerVector& Fvec,
                            const double& Jss1,
                            const double& Js1s, 
                            const double& P0,
                            const int&  minimum_n,
                            const double& start_point=0.0,
                            const double& end_point=0.0);

#endif