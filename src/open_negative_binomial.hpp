#ifndef open_negative_binomial_header
#define open_negative_binomial_header

using namespace Rcpp;

double prob_approx( const double& count,
                    const double& split_U,
                    const double& alpha,
                    const double& mu,
                    const int& minimum_n,
                    const double& tau0,
                    const double& tau1,
                    const double& error);

#endif
