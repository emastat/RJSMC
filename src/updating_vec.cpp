#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//' @title augment_vecZQFV
//' Function for augment/updating the state vectors (Zvec, Vvec, Fvec) as resulting from accepting a forward step
//' @param vec (IntegerVector&) vector to be augmented
//' @param left_value (const int) state value of the new left segment
//' @param right_value (const int) state value of the new right segment
//' @param IN (const int) index of the the  segment splitted(range is (1,S) )
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void augment_vecZQFV(IntegerVector& vec,
                     const int left_value,
                     const int right_value,
                     const int IN,
                     const int S){

  // store the new left state

    vec[IN-1] = left_value ;

 // move state value from IN+1 to S one cell to the right

    for(int j= S ; j>IN ;j-- ){


      vec[j] = vec[j-1] ;

    }

    vec[IN] = right_value ;


}

//' @title augment_vecM
//' Function for augment/updating the state vector Mvec as resulting from accepting a forward step
//' @param vec (NumericVector&) vector to be augmented
//' @param left_value (const double) state value of the new left segment
//' @param right_value (const double) state value of the new right segment
//' @param IN (const int) index of the  segment splitted (range is (1,S) )
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void augment_vecM(NumericVector& vec,
                  const double left_value,
                  const double right_value,
                  const int IN,
                  const int S){

  // store the new left state

  vec[IN-1] = left_value ;

  // move state value from IN+1 to S one cell to the right

  for(int j= S ; j>IN ;j-- ){


    vec[j] = vec[j-1] ;

  }

  vec[IN] = right_value ;

}


//' @title reduce_vecZQFV
//' Function for reducing/updating the state vectors (Zvec, Vvec, Fvec) as resulting from accepting a backward step
//' @param vec (IntegerVector&) vector to be augmented
//' @param  (const int) state value of the new segment
//' @param BR (const int) index of the removed breakpoint (range (1,S-1) )
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void reduce_vecZQFV(IntegerVector& vec,
                    const int value,
                    const int BR,
                    const int S){

  // store the new left state

  vec[BR-1] = value ;

  // move state value from IN+1 to S one cell to the right

  for(int j= BR+1 ; j<S ;j++ ){


    vec[j-1] = vec[j] ;

  }

  // set the ¨ of the last segment to 0

  vec[S-1] = 0 ;

}


//' @title reduce_vecM
//' Function for reducing/updating the state vector Mvec as resulting from accepting a backward step
//' @param vec (NumericVector&) vector to be augmented
//' @param  (const double) state value of the new segment
//' @param BR (const int) index of the removed breakpoint (range (1,S-1) )
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void reduce_vecM(NumericVector& vec,
                  const double value,
                  const int BR,
                  const int S){

  // store the new left state

  vec[BR-1] = value ;

  // move state value from IN+1 to S one cell to the right

  for(int j= BR+1 ; j<S ;j++ ){


    vec[j-1] = vec[j] ;

  }

  // set the ¨ of the last segment to 0

    vec[S-1] = 0.0 ;

}



//' @title augment_Bvec
//' Function for augmentig the breakpoint vector Bvec as result from accepting a forward step
//' @param vec (NumericVector&) vector with breakpoint
//' @param T_p (const double) vector to be inserted
//' @param IN (const int) index of the segment splitted
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void augment_Bvec(NumericVector& Bvec,
                  const double T_p,
                  const int IN,
                  const int S){


  // move 1 position to the right all the breakpoint from IN onwards

  for(int j=S+1 ; j>IN ; j--){

    Bvec[j] = Bvec[j-1] ;

  }

  Bvec[IN] = T_p ;

}



//' @title reduce_Bvec
//' Function for reducing the breakpoint vector Bvec as result from accepting a backward step
//' @param vec (NumericVector&) vector with breakpoint
//' @param BR (const int) index of the breakpoint to be removed
//' @param S (const int) number of current segments
//' @return void no output is returned because the input (vec) is passed by reference and hence modified
//' @export
// [[Rcpp::export]]
void reduce_Bvec(NumericVector& Bvec,
                 const int BR,
                 const int S){


  // shift all the breakpoints from index BR+1 onwards 1 position to the left

   for(int j=BR+1 ; j<S+2;j++){

     Bvec[j-1] =Bvec[j] ;

   }

   //Bvec[S] = 0.0 ;
}











