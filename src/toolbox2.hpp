#ifndef toolbox2_cpp_header
#define toolbox2_cpp_header

using namespace Rcpp;

NumericVector stabilize(const NumericVector weight_vec);

IntegerVector table_cpp(IntegerVector X, int n,bool include_zero);

#endif