#ifndef tool_functions_cpp_header
#define tool_functions_cpp_header

using namespace Rcpp;

int r_rcatlp(int n, NumericVector log_prob);

double prod(NumericVector x_vec);

NumericVector make_table( IntegerVector X, int K,bool inc_zero=false);


NumericVector reptest2(NumericVector x, NumericVector y) ;


NumericVector r_rdirichlet(int n, NumericVector alpha);


double r_dtpois(double x, double lambda);

double r_rtgamma(double shape, double rate, double truncation);

IntegerVector cpp_findInterval(NumericVector x, NumericVector breaks);

int factorial_cpp(int n);

double lnfactorial( int a);

CharacterVector cpp_paste(CharacterVector lhs, CharacterVector rhs);

IntegerVector which_true(LogicalVector y );

IntegerVector cpp_order(NumericVector x) ;

NumericVector cpp_rep(NumericVector x, NumericVector y);

NumericVector smart_round( NumericVector X,  int digi );

NumericVector seq_cpp(double start, double end, double gap);

int not_null(const List& x);

List which_not_null(const List& X);

double rcpp_sum(const NumericVector& X);

#endif
