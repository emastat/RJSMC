#ifndef RJMCMC_SMC_header
#define RJMCMC_SMC_header

using namespace Rcpp;


void     RJMCMC_SMC(const NumericVector& T_seg,
                    const IntegerVector& Y_seg,
                    const int& U,
                    const int& K,
                    const int& W,
                    const double& start_point,
                    const double& end_point,
                    const double& t_star,
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
                    const double& P0,
                    const int& minimum_n,
                    const double& Jss1,
                    const double& Js1s,
                    const int& Smax,
                    const int& n_ite,
                    const int& burn_in,
                    const int& thinning,
                    const int& n_particle,
                    const IntegerVector& particle_index_vec,
                    IntegerVector& V_last_complete,
                    NumericVector& B_last,
                    List& container_B,
                    List& container_V,
                    List& container_Z,
                    List& container_Q,
                    List& container_F,
                    IntegerVector& Svec,
                    NumericVector& weight_vec,
                    const bool& empty_seg);

#endif
