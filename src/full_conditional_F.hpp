#ifndef full_conditional_F_header
#define full_conditional_F_header

using namespace Rcpp;

List full_conditional_F(const int& F,
                        const double& L,
                        const NumericVector& key0vec,
                        const NumericVector& eta0vec,
                        const NumericVector& probvec_F,
                        const bool sample_F=false);

#endif