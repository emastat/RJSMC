#ifndef weigths_correction_header
#define weigths_correction_header


using namespace Rcpp;

double weigths_correction(const NumericVector& T_seg,
                          const List& break_sample_output,
                          const int& M,
                          const int& K,
                          const NumericVector& keyvec,
                          const NumericVector& etavec,
                          const NumericVector& key0vec,
                          const NumericVector& eta0vec,
                          const NumericVector& probvec_Z,
                          const NumericVector& probvec_F,
                          const double& P0);

#endif
