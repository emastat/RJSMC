#ifndef updating_vec_header
#define updating_vec_header


using namespace Rcpp;

void augment_vecZQFV(IntegerVector& vec,
                     const int left_value,
                     const int right_value,
                     const int IN,
                     const int S);

void augment_vecM(NumericVector& vec,
                  const double left_value,
                  const double right_value,
                  const int IN,
                  const int S);


void reduce_vecZQFV(IntegerVector& vec,
                    const int value,
                    const int BR,
                    const int S);

void reduce_vecM(NumericVector& vec,
                 const double value,
                 const int BR,
                 const int S);


void augment_Bvec(NumericVector& Bvec,
                  const double T_p, 
                  const int IN,
                  const int S);


void reduce_Bvec(NumericVector& Bvec,
                 const int BR,
                 const int S);


#endif