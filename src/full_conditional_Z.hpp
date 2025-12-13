#ifndef full_conditional_Z_header
#define full_conditional_Z_header

using namespace Rcpp;

List full_conditional_Z(const int& Z,
                        const int& K,
                        const double& L,
                        const NumericVector& keyvec,
                        const NumericVector& etavec,
                        const NumericVector& probvec_Z,
                        const bool& sample_Z=false);


#endif

