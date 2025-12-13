#include <Rcpp.h>

using namespace Rcpp;

//' @title systematic_resampling:  function for performing systematic resampling
//'
//' @param weight_vec vector with the current normalized weight for each particle
//' @param container_B List for with Bvec generated inside the current Update Interval, for each particle
//' @param container_V List for with Vvec generated inside the current Update Interval, for each particle
//' @param container_Z  List for with Zvec generated inside the current Update Interval, for each particle
//' @param container_Q  List for with Qvec generated inside the current Update Interval, for each particle
//' @param container_F  List for with Fvec generated inside the current Update Interval, for each particle
//' @param Svec vector with number of segment in each particle
//' @param V_last_complete Vector with the V state for the last complete segment of each particle
//' @param B_last (NumericVector&) Vector with the last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration
//' @param n_particle number of particles created
//' @return The function does not explicitly return any object: However, by reference 9 element are updated:
//' @export
// [[Rcpp::export]]
List resampling_func(const NumericVector& weight_vec,
                     const List& container_B,
                     const List& container_V,
                     const List& container_Z,
                     const List& container_Q,
                     const List& container_F,
                     const IntegerVector& Svec,
                     const IntegerVector& V_last_complete,
                     const NumericVector& B_last,
                           int n_particle){


 // sample first cut off of the systematic resampling
 double u1 = R::runif(0,1.0/n_particle) ;


 // import the R function seq
 Function seq_r("seq");

 // generate vector with cut offs
 NumericVector u_vec = seq_r(Named("from")=u1, Named("to")=1.0, Named("by")=1.0/n_particle);

 // cumulative density for the weight vector
 NumericVector sum_weight_vec = cumsum(weight_vec) ;

 // initialize required values
 int current_u = 0 ;
 int current_particle = 0 ;
 double previous_sum_w = 0.0 ;
 IntegerVector new_index_vec(n_particle);

 double current_sum_w = sum_weight_vec[current_particle] ;

 while(current_u < n_particle){


   if( (previous_sum_w<u_vec[current_u]) & (current_sum_w>=u_vec[current_u] )){

     // we have the particle occupying the current slot

     new_index_vec[current_u] = current_particle;

     //move to the next slot available

     current_u = current_u + 1 ;

    }else{

     // current particle cannot fill the slot
     // move to next particle

     previous_sum_w = current_sum_w ;
     current_particle = current_particle + 1 ;
     current_sum_w = sum_weight_vec[current_particle] ;

   }
  }


 // update all the particles using nee container (prefix sr=systematic resampling)

 int new_idx;
 List sr_container_B(n_particle); // List for storing the Bvec generated inside the current Update Interval
 List sr_container_V(n_particle); // List for storing the Vvec generated inside the current Update Interval
 List sr_container_Z(n_particle); // List for storing the Zvec generated inside the current Update Interval
 List sr_container_Q(n_particle); // List for storing the Qvec generated inside the current Update Interval
 List sr_container_F(n_particle); // List for storing the Fvec generated inside the current Update Interval
 IntegerVector sr_Svec(n_particle) ; // number of segments currently in each particle

 // last simulated changepoint WITHIN the Update Interval, at the previous SMC iteration
 NumericVector sr_B_last = rep(0.0,n_particle) ;
 IntegerVector sr_V_last_complete = rep(1,n_particle) ; // V state for the last complete segment of each particle

 for(int i=0; i<n_particle; i++){

   new_idx= new_index_vec[i] ;

   sr_container_B[i] = clone(as<NumericVector>(container_B[new_idx]));
   sr_container_V[i] = clone(as<NumericVector>(container_V[new_idx]));
   sr_container_Z[i] = clone(as<NumericVector>(container_Z[new_idx]));
   sr_container_Q[i] = clone(as<NumericVector>(container_Q[new_idx]));
   sr_container_F[i] = clone(as<NumericVector>(container_F[new_idx]));
   sr_Svec[i] = Svec[new_idx];
   sr_B_last[i] = B_last[new_idx];
   sr_V_last_complete[i] = V_last_complete[new_idx] ;

 }

 return (List::create(
     Named("sr_container_B")=sr_container_B,
     Named("sr_container_V")=sr_container_V,
     Named("sr_container_Z")=sr_container_Z,
     Named("sr_container_Q")=sr_container_Q,
     Named("sr_container_F")=sr_container_F,
     Named("sr_Svec")= sr_Svec,
     Named("sr_V_last_complete")= sr_V_last_complete,
     Named("sr_B_last")= sr_B_last
 ));
}
