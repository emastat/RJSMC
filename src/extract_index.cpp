#include <Rcpp.h>


using namespace Rcpp;

//' @title extract_index
//' Function for extracting the first and last index elements of an ordered vector
//' sastisfying a constraint specified by a lower and upper bound
//' @param  x (const NumericVector&) the original vector
//' @param  lower (const double), the lower bound (included in the subsetting)
//' @param  upper (const double), the lower bound
//' @return 3-element vector with:
//' \describe{
//' \item{element [0]}{length of the subsetted vector}
//' \item{element [1]}{index of the first element of the subsetted vector }
//' \item{element [2]}{index of the last element of the subsetted vector. If the last index is "-1" then no elements sastisfy the constraint}
//' }
// [[Rcpp::export]]
NumericVector extract_index(const NumericVector& x,
                            const double& lower,
                            const double& upper) {

  NumericVector out=(3) ;

   int i = 0;

   int x_length = x.size() ;

  while((x[i] <= lower) & (i < x_length)) i = i + 1;

  int j = i;

  out[1] = i ; //store the left index

 if(j<x_length){

    while((x[j] < upper)  & (j < x_length) ) j = j + 1;

  j=j-1 ;

  if(j>=i){

    out[0] = j-i +1 ; //store the length of the selected vector

    out[2] = j ; //store the right index

  }


 }else{ j =-1 ; out[2] = j ; }



  return  out ;
}




