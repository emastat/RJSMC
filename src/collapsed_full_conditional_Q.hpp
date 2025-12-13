#ifndef collapsed_full_conditional_Q_header
#define collapsed_full_conditional_Q_header

using namespace Rcpp;



List collapsed_full_conditional_Q(const int& Q,
                                  const int& W,
                                  const double& count,
                                  const double& LB,
                                  const double& UB,
                                  const NumericVector& alphavec,
                                  const NumericVector& muvec,
                                  const NumericVector& probvec_Q,
                                  const int& minimum_n,
                                  const bool open_segment=false,
                                  double end_point= 0.0,
                                  const bool sample_Q=false);

#endif