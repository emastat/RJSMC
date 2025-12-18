#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

int r_rcatlp(int n, NumericVector log_prob){

  Rcpp::Environment extraDistr("package:extraDistr");

  //calling rcatlp
  Function f = extraDistr["rcatlp"];
  int  value = as<int>(f(_["n" ]=n, _["log_prob"]=log_prob));
  value+=1 ;

  return value;
}


double prod(NumericVector x_vec){

  double prod_value=1;

  for(int i=0;i<x_vec.length();i++){

   prod_value *= x_vec[i] ;
  }

  return prod_value;
}

//' @title Function for creating ax frequency table.
//' @param X NumericVector. Vector whose element are used to create the frequency table
//' @param K Integer. The number of classes -max(X) standard choice-
//' @param inc_zero Boolean. inc_zero==F (default) if the 0's must not be taken into account in the counting.
//' @return The frequency vector (NumericVector)
// [[Rcpp::export]]
NumericVector make_table( IntegerVector X, int K,bool inc_zero=false){

  NumericVector freque(K+1);


  LogicalVector to_sum ;

  for(int i=0; i<K+1;i++){

    to_sum= X==i  ;
    freque[i] = Rcpp::sum(to_sum) ;

  }

  if(inc_zero==false){

    freque.erase(0);
  }

  return freque;
}


NumericVector reptest2(NumericVector x, NumericVector y) {
  int n = y.size();
  NumericVector myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
    ind += p;
  }
  return myvector;
}



NumericVector r_rdirichlet(int n, NumericVector alpha){

  Rcpp::Environment extraDistr("package:extraDistr");
  //calling dtpois
  Function f = extraDistr["rdirichlet"];
  NumericVector value = as<NumericVector>(f(_["n" ]=n,_["alpha"]=alpha) );

  return value;
}


double r_dtpois(double x, double lambda){

  Rcpp::Environment extraDistr("package:extraDistr");
  //calling dtpois
  Function f = extraDistr["dtpois"];
  double  value = as<double>(f(_["x" ]=x,_["lambda"]=lambda,_["a"]=0,_["log"]=true));

  return value;
}


 double r_rtgamma(double shape, double rate, double truncation){

   Rcpp::Environment truncdist("package:truncdist");

   //calling rtrunc

   Function f = truncdist["rtrunc"];

   double value = as<double>( f(_["n"]=1, _["spec"]="gamma", _["a"]=truncation, _["b"]=R_PosInf, _["shape"]=shape, _["rate"]=rate) );

  return value ;

 }


// [[Rcpp::plugins(cpp11)]]

//' @title Function to assign each element of a vector "x" to one of N  unique intervals defined through the
//' vector of breakpoint B; \code{B[j]<=x<B[j+1]} j=1...N+1. The out is "j"
//' @param x NumericVector. Vector  whose element need to ble allocated
//' @param breaks NumericVector. Vector containing the N+1 breakpoints
//' @return IntegerVector having the same size as x. Each element contains the index referring to the interval \code{x[i]} belongs to
//' @export
// [[Rcpp::export]]
IntegerVector cpp_findInterval(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());

  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;

  int max_index = breaks.size() - 1;

  for(it = x.begin(), out_it = out.begin(); it != x.end();
  ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    int index = std::distance(breaks.begin(), pos);
    *out_it = (index > max_index) ? max_index : index;
  }

  return out;
}

// C++ program to find factorial of given number

int factorial_cpp(int n) {
  int res = 1, i;
  for (i = 2; i <= n; i++)
    res *= i;
  return res;
}

double lnfactorial( int a){

  int y;
  double z;
  if (a == 1)
    return 0;
  else
  {
    z = 0;

    for (y = 2; y<=a; y++ )

      z = log(y)+z;


    return z;
  }
}

CharacterVector cpp_paste(CharacterVector lhs, CharacterVector rhs){
  using proxy_t = internal::string_proxy<STRSXP>;

  std::vector<std::string> res(lhs.begin(), lhs.end());
  std::transform(res.begin(), res.end(), rhs.begin(), res.begin(),
                 [&](const std::string& x, const proxy_t& y) {
                   return x + y;
                 }
  );

  return wrap(res);
}

IntegerVector which_true(LogicalVector y ){

  using namespace Rcpp;

  // [[Rcpp::plugins(cpp11)]]

  IntegerVector index(y.size()) ;

  int k=0 ;

  for(int i=0; i<y.length();i++){


    if(y[i]==true){   index[k]=i ; k=k+1  ;}

  }

  index.erase(k,y.size()) ;

  return index;

}

IntegerVector cpp_order(NumericVector x){

  NumericVector sorted = clone(x).sort();

  IntegerVector value = match(sorted, x) -1;
  return value;
}

NumericVector cpp_rep(NumericVector x, NumericVector y){
  int n = y.size();
  NumericVector myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(myvector.begin()+ind, myvector.begin()+ind+p, x[i]);
    ind += p;
  }
  return myvector;
}

NumericVector smart_round( NumericVector X,  int digi ){


  double up = std::pow(10,digi);

  NumericVector x = X * up ;

  NumericVector y = Rcpp::floor(x);


  IntegerVector ordered_xy = cpp_order(x-y) ;

  double sum_x = Rcpp::sum(x);

  double sum_y = Rcpp::sum(y);

  double  round_sum = std::round( sum_x - sum_y ) ;


  IntegerVector indices = tail( ordered_xy,round_sum ) ;

  NumericVector  add_one(x.size()) ;

  add_one[indices] = 1.0 ;

  y = y + add_one ;


  NumericVector value = y / up ;


  return value ;

}

NumericVector seq_cpp(double start, double end, double gap){

  int len = (end - start )/gap + 1  ;

  NumericVector out(len ) ;

  out[0] = start ;
  int i = 1  ;

  while(start<=end-gap){

    out[i] = start + gap ;

    start = start + gap ;

    i = i + 1 ;

  }

  return out ;

}

int not_null(const List& x){

  bool out = x[0]!=R_NilValue ;

  return out ;

}

List which_not_null(const List& X){


  IntegerVector Y = sapply( X, not_null );

  List W = X[Y==1] ;
  return W ;

}

double rcpp_sum(const NumericVector& X){

  double out = sum(X);

  return out ;

}











