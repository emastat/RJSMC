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
//' \item{element \code{[0]}}{length of the subsetted vector}
//' \item{element \code{[1]}}{index of the first element of the subsetted vector }
//' \item{element \code{[2]}}{index of the last element of the subsetted vector. If the last index is "-1" then no elements sastisfy the constraint}
//' }
// [[Rcpp::export]]
NumericVector extract_index(const NumericVector& x,
                            const double& lower,
                            const double& upper) {

  NumericVector out=(3) ;

   int i = 0;
   int x_length = x.size() ;

  while((i < x_length) && (x[i] <= lower)) {
    i = i + 1;
  }

  // If all elements are <= lower, no valid start index exists
  if(i >= x_length){
    out[0] = 0;  // count = 0
    out[1] = -1; // no valid start index
    out[2] = -1; // no valid end index
    return out;
  }

  int j = i;
  out[1] = i ;

 if(j<x_length){

    while((j < x_length) && (x[j] < upper)) {
      j = j + 1;
    }

    j = j - 1 ;
    
    if(j < 0){
      j = -1;
      out[2] = j;
      return out;
    }
    if(j >= x_length){
      j = x_length - 1;
    }

    if(j>=i){
      out[0] = j-i +1 ;
      out[2] = j ;
    }


 }else{ j =-1 ; out[2] = j ; }



  return  out ;
}




