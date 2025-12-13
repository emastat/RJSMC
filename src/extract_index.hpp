#ifndef extract_index_header
#define extract_index_header

using namespace Rcpp;



NumericVector extract_index(const NumericVector& x, 
                            const double& lower, 
                            const double& upper);

#endif