#ifndef initial_breakpoints_header
#define initial_breakpoints_header

using namespace Rcpp;



NumericVector initial_breakpoints(  const NumericVector& T_seg_orig,
                                    const double& t_star,
                                    const double& start_point,
                                    const double& end_point,
                                    const int& minimum_n,
                                    const double& max_range,
                                    const double& min_seg_length);


#endif
