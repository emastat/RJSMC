#include <Rcpp.h>
#include "toolbox1.hpp"

using namespace Rcpp;

//' @title initial_breakpoints
//' Function to initialize the changepoint vector during the first iteration of the RJMCMC
//' @param T_seg_orig Vector with time stamps within the update interval
//' @param t_star Estimate of the latest changepoint simulated when performing the RJMCMC during the previous iteration of the SMC (precedent update interval)
//' @param start_point Left bound of the current update interval
//' @param end_point Right bound of the current update interval
//' @param minimum_n minimum number of observation allowed in a non-empty segment
//' @param max_range (const double&) time interval in which the last changepoint is allowed to be place after the end point (end_point;end_point+1)
//' @param min_seg_length (const double) minimum length allowed for a segment to be valid
//' @return NumericVector with the generated changepoints. Starting with t_star, ending with B_k > UB
//' @export
// [[Rcpp::export]]
NumericVector initial_breakpoints(  const NumericVector& T_seg_orig,
                                    const double& t_star,
                                    const double& start_point,
                                    const double& end_point,
                                    const int& minimum_n,
                                    const double& max_range,
                                    const double& min_seg_length){

  //remove time stamps smaller than t_tar
  NumericVector T_seg = clone(T_seg_orig) ;

  T_seg = T_seg[T_seg>t_star] ;

  // generate and order the breakpoints

  if(T_seg.size()>=minimum_n){

    int T_seg_length = T_seg.length();

    // T_segplus: vector with observed time-stamps augmented with t_star (in 0)
    NumericVector T_segplus(T_seg_length+1) ;

    T_segplus[0] = t_star ;

    T_segplus[ Range(1,(T_segplus.size() - 1 )) ] = T_seg  ;

    // int_vec: vector with times between each time-stamps

    NumericVector  int_vec = diff(T_segplus) ;

    //initialize group vec - it will contain the group index each time-stamp will be assigned to

    IntegerVector group_vec = rep(0, T_seg_length) ;


    //observations in (t_star,start_point) must belong to group "1"

    int i = 1 ;

    group_vec[T_seg<start_point] =  i ;

    int starting_point = sum(T_seg<start_point) ;


    // assign sequentially each observation to the current or the next segment

    for(int j=starting_point ; j<T_seg_length;j++){

      //select the "time diff btw two time-stamps" assigned to group "i"

      NumericVector int_vec_group_i =  int_vec[group_vec==i] ;

      //do the mean of these selected time spaces

      double ave =  0.0 ;

      if(int_vec_group_i.size()>0){

        ave =  mean(int_vec_group_i) ;

      }else{ ave = 0.0 ; }

      //make sure that each group has at least minimum_n messages &
      //the segment is long at least "min_seg_length"

      if(sum(group_vec==i) < minimum_n){

        group_vec[j] = i ;  //assign obs "j" to group "i"

      }else if(sum(int_vec_group_i) < min_seg_length){

        group_vec[j] = i ;  //assign obs "j" to group "i"

      }else if(ave < 0.0045){

        group_vec[j] = i ;

      }else if(int_vec[j]<=ave*1.5 ){

        group_vec[j] = i ;

      }else{

        i = i + 1 ;

        group_vec[j] = i ;

      }

    }


    // diff_group: vector of the same length of T_seg where each element is 0
    //unless it reflects the position of the last time-stamp of each segment (in which case the value is 1)

    group_vec.push_front(1) ;

    IntegerVector diff_group = diff(group_vec) ;


    double cost ;

    // index_group: vector with the position of the last time-stamp of each segment
    IntegerVector index_group = which_true(diff_group==1);

    NumericVector T_point(2 * index_group.size() ) ;

    int in0 = 0 ; int in1 = 1 ;

    for(int i : index_group){

      cost    = int_vec[i]*0.02 ;

      T_point[in0]  = T_seg[i-1] + cost ;

      T_point[in1] = T_seg[i] - cost ;

      in0 = in1 + 1 ; in1 = in1 + 2 ;

    }

    NumericVector Bvec1(T_point.size() + 1,-100) ;


    //update the last breakpoint: always after the end of the Update interval
    Bvec1[Bvec1.size()-1] = R::runif(end_point,end_point+max_range) ;

    if(T_point.size()>0){

      Bvec1[Range(0, Bvec1.size()-2  )] = T_point ;

    }


    //remove breakpoints smaller than start_point
    NumericVector Bvec1_final = Bvec1[Bvec1>start_point] ;

    //add t_star as first breakpoint
    Bvec1_final.push_front(t_star) ;
    
    // remove breakpoints creating segments smaller than 5 seconds

    NumericVector seg_length = diff(Bvec1_final) ;


    // check that no closed segment has less than minimum_n observations

    IntegerVector  obs_position = cpp_findInterval(T_seg,Bvec1_final);

    IntegerVector segment_table = table(obs_position) ;

    int previous_seg = 1;

    int current_seg = 1;

    for( int i=0 ; i< (segment_table.size() - 2) ; i++){

      previous_seg = current_seg ;

      int count_value = segment_table[i] ;

      if(count_value >= minimum_n){

        current_seg = 1 ;

      }else if((count_value>0) & (count_value< minimum_n)){

        Rcout << "obs_position: " << obs_position << "\n" ;
        Rcout << "Bvec1_final: " << Bvec1_final << "\n" ;
        Rcout << "T_seg: " << T_seg << "\n" ;
        Rcout << "t_star: " << t_star << "\n" ;
        Rcout << "start_point: " << start_point << "\n" ;

        stop("number of obs in a closed segment must be 0 or >minimum_n: initial_breakpoint");

      }else{

        if(previous_seg == 0){stop("both the previous and the current segment are null") ;}

        current_seg = 0 ;

      }

    }

    bool check_order = std::is_sorted(Bvec1_final.begin(),Bvec1_final.end()) ;

    if(check_order == false){

      Rcout << "Initial_breakpoint: this is Bvec:  " << Bvec1_final << "\n"  ;

      stop("Initial_breakpoint: the Breapoint vector is not sorted");


      }


    return Bvec1_final ;

  }else if((T_seg.size()>0)  && (T_seg.size() < minimum_n)){

    //if less than minimum_n (but more than 0), create only 1 segment

    NumericVector Bvec1_final(2) ;

    Bvec1_final[0] = t_star ;
    Bvec1_final[1] =R::runif(end_point,end_point+ max_range) ;


    return Bvec1_final ;

  }else{

    stop("0 observation has been found when creating the initial breakpoints vector: it's not possible");

  }


}

