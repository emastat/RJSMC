#include <Rcpp.h>
#include "toolbox2.hpp"
#include "extract_index.hpp"
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
//' @title full_conditional_V
//' Function for sampling the V  state for a  new non-empty segment from its full conditional and/or evaluate it at the sampled point
//' @param V (const int&) Value of the current V if present; if not set it 0.
//' @param U (const int&) Number of levels for the V state
//' @param num_logs (const int&) Total number of unique messages that can be observed
//' @param T_seg (const NumericVector&) Vector with time stamps within the update interval
//' @param Y_seg (const IntegerVector&) Vector with the messages within the update interval
//' @param LB (const double&) lower bound of the selected segment
//' @param UB (const double&) upper bound of the selected segment
//' @param lambdamat (const NumericMatrix&) probability weights ruling the message composition of the segments for the different levels of the V state
//' @param probvec_V (const NumericVector&) Probability mass for the V state distribution
//' @param V_left  V value for the previous segment
//' @param open_segment (const bool, default is false) if false  the segment is considered totally observed. If true the semgent is considered partially observed; this affects the likelihood computation
//' @param end_point (double, default is 0.0) to be used for partially observed segments; declares the point where the stops to be observed
//' @param P0 (const double&) Probability to observe an empty segment after a non-empty one
//' @param sample_V ( const bool=false) if true a new V is sampled and the full conditional is evaluated at the new value. If false (by default) the full condititional is evaluated at the input V
//' @return A list containting 2 elements:
//' \describe{
//' \item{V_new}{int, the sample Vstate for the selected segment or the one inputed}
//' \item{eval_densV}{double, value of the log full conditional for either the new sampled V or the current Q inputed ( depending on sample_V)}
//' }
//' @name full_conditional_V
// [[Rcpp::export]]
List full_conditional_V(const int& V,
                        const int& U,
                        const int& num_logs,
                        const NumericVector& T_seg,
                        const IntegerVector& Y_seg,
                        const double& LB,
                        const double& UB,
                        const NumericMatrix& lambdamat,
                        const NumericVector& probvec_V,
                        int V_left,
                        const double& P0,
                        const double& end_point=0.0,
                        const bool open_segment=false,
                        const bool sample_V=false){


  // initialize V and other the required quantities

  int V_new=0 ;

  double eval_densV, weight, prior_V ;

  NumericVector log_densV(U), densV(U);

  NumericVector prod_dj_logLL(num_logs) ;


  //initialize out_extract: vector with indexes of the observations contained in the segment
  NumericVector out_extract(3) ;

  // compute log_vector: vector with messages within the selected segment

  if(open_segment==false){

    // observations found in the whole segment
    out_extract =clone( extract_index(T_seg,LB,UB)) ;

  }else{

    // observations found in the partial segment
    out_extract =clone( extract_index(T_seg,LB,end_point)) ;
  }

  IntegerVector log_vector;
  if(out_extract[2] < 0){
    log_vector = IntegerVector(0) ;
  } else {
    int end_idx = (int)out_extract[2];
    int start_idx = (int)out_extract[1];
    
    if(end_idx >= Y_seg.size()){
      end_idx = Y_seg.size() - 1;
    }
    if(start_idx < 0){
      start_idx = 0;
    }
    if(start_idx > end_idx){
      log_vector = IntegerVector(0) ;
    } else {
      log_vector = Y_seg[Range(start_idx, end_idx)] ;
    }
  }



  // frequency table for the messages
  IntegerVector dj_vec = table_cpp(log_vector,num_logs,false) ;

  // Extract count of observations for validation check
  double count = out_extract[0] ;
  int start_idx_extract = (int)out_extract[1];
  int end_idx_extract = (int)out_extract[2];
  
  // ENHANCED BUG DETECTION: Directly count observations in the segment
  int actual_count = 0;
  if(open_segment == true){
    // For open segments, count observations in [LB, end_point)
    for(int k = 0; k < T_seg.size(); k++){
      if(T_seg[k] >= LB && T_seg[k] < end_point){
        actual_count++;
      }
    }
  } else {
    // For closed segments, count observations in [LB, UB)
    for(int k = 0; k < T_seg.size(); k++){
      if(T_seg[k] >= LB && T_seg[k] < UB){
        actual_count++;
      }
    }
  }
  
  // DEBUG: If counts don't match, print detailed information about extract_index results
  if(open_segment == true && count != actual_count){
    Rcout << "\n=== DEBUG extract_index mismatch ===\n";
    Rcout << "  extracted_count=" << count << ", actual_count=" << actual_count << "\n";
    Rcout << "  LB=" << LB << ", end_point=" << end_point << ", UB=" << UB << "\n";
    Rcout << "  extract_index returned: count=" << out_extract[0] << ", start_idx=" << out_extract[1] << ", end_idx=" << out_extract[2] << "\n";
    Rcout << "  T_seg size=" << T_seg.size() << "\n";
    
    // Manually trace through extract_index logic
    Rcout << "\n  Manual trace of extract_index logic:\n";
    int i_trace = 0;
    int x_length_trace = T_seg.size();
    Rcout << "    Step 1: Finding first index where T_seg[i] > " << LB << "\n";
    while((i_trace < x_length_trace) && (T_seg[i_trace] <= LB)) {
      Rcout << "      i=" << i_trace << ", T_seg[" << i_trace << "]=" << T_seg[i_trace] << " <= " << LB << ", continue\n";
      i_trace = i_trace + 1;
    }
    if(i_trace >= x_length_trace){
      Rcout << "    i reached end, capping to " << (x_length_trace - 1) << "\n";
      i_trace = x_length_trace - 1;
    }
    Rcout << "    Final i=" << i_trace << ", T_seg[" << i_trace << "]=" << T_seg[i_trace] << "\n";
    
    int j_trace = i_trace;
    Rcout << "    Step 2: Finding first index where T_seg[j] >= " << end_point << "\n";
    if(j_trace < x_length_trace){
      while((j_trace < x_length_trace) && (T_seg[j_trace] < end_point)) {
        Rcout << "      j=" << j_trace << ", T_seg[" << j_trace << "]=" << T_seg[j_trace] << " < " << end_point << ", continue\n";
        j_trace = j_trace + 1;
      }
      Rcout << "    j after while loop=" << j_trace << "\n";
      j_trace = j_trace - 1;
      Rcout << "    j after decrement=" << j_trace << "\n";
      
      if(j_trace < 0){
        Rcout << "    j < 0, setting to -1\n";
      } else if(j_trace >= x_length_trace){
        Rcout << "    j >= x_length, capping to " << (x_length_trace - 1) << "\n";
        j_trace = x_length_trace - 1;
      }
      
      if(j_trace >= i_trace){
        int computed_count = j_trace - i_trace + 1;
        Rcout << "    j >= i, computed count = " << j_trace << " - " << i_trace << " + 1 = " << computed_count << "\n";
        if(j_trace >= 0 && j_trace < x_length_trace){
          Rcout << "    T_seg[" << j_trace << "]=" << T_seg[j_trace] << " (should be < " << end_point << ")\n";
        }
        if(i_trace >= 0 && i_trace < x_length_trace){
          Rcout << "    T_seg[" << i_trace << "]=" << T_seg[i_trace] << " (should be > " << LB << ")\n";
        }
      } else {
        Rcout << "    j < i, no count computed\n";
      }
    }
    
    // Print observations in the range
    Rcout << "\n  Observations in range [" << LB << ", " << end_point << "):\n";
    int obs_in_range = 0;
    for(int k = 0; k < T_seg.size(); k++){
      if(T_seg[k] >= LB && T_seg[k] < end_point){
        Rcout << "    T_seg[" << k << "]=" << T_seg[k] << "\n";
        obs_in_range++;
      }
    }
    Rcout << "  Total observations in range: " << obs_in_range << "\n";
    
    // Check what extract_index actually found
    if(end_idx_extract >= 0 && end_idx_extract < T_seg.size()){
      Rcout << "\n  extract_index claims observation at index " << end_idx_extract << ": T_seg[" << end_idx_extract << "]=" << T_seg[end_idx_extract];
      if(T_seg[end_idx_extract] >= LB && T_seg[end_idx_extract] < end_point){
        Rcout << " (VALID - in range)\n";
      } else {
        Rcout << " (INVALID - NOT in range!)\n";
      }
      if(start_idx_extract >= 0 && start_idx_extract < T_seg.size()){
        Rcout << "  extract_index start index " << start_idx_extract << ": T_seg[" << start_idx_extract << "]=" << T_seg[start_idx_extract];
        if(T_seg[start_idx_extract] >= LB && T_seg[start_idx_extract] < end_point){
          Rcout << " (VALID - in range)\n";
        } else {
          Rcout << " (INVALID - NOT in range!)\n";
        }
      }
    }
    Rcout << "=== END DEBUG ===\n\n";
  }

  // compute probability weights for the full conditional

  for( int i=0; i<U ; i++){

    for(int j=0 ; j<num_logs ; j++){

      prod_dj_logLL[j] =  std::log( lambdamat(i,j) ) * (double)dj_vec[j] ;

    }

    weight = sum(prod_dj_logLL) ;

    if( V_left==0){

      prior_V = probvec_V[i] ;

    }else{

      prior_V =(1-P0) * probvec_V[i] ;

    }

    log_densV[i] = std::log(prior_V) + weight  ;

  }


  // employ the log-sum-exp trick

  densV = stabilize(log_densV) ; // final probability vector


  if(sample_V==true){

    // a new Q must be sampled - log full conditional evaluated in this point

    // SAMPLE V
    V_new = as<int>( Rcpp::sample( U, 1, false, densV, true ) ) ;

    // logarithm of the evaluated the density in the sampled point
    eval_densV = std::log(densV[V_new - 1 ]) ;

    // ENHANCED BUG DETECTION CHECK: Open segment with zero observations must have V = 0
    // Check both the extracted count and the actual count
    if(open_segment == true && (count == 0.0 || actual_count == 0) && V_new != 0){
      Rcpp::stop("BUG DETECTED in full_conditional_V: Open segment with zero observations has non-zero V state. V=%d, extracted_count=%.0f, actual_count=%d, start_idx=%d, end_idx=%d, LB=%.6f, UB=%.6f, end_point=%.6f, T_seg_size=%d", 
                 V_new, count, actual_count, start_idx_extract, end_idx_extract, LB, UB, end_point, T_seg.size());
    }

    // return a list with the new sample Q and the value of the full conditional in this value

    return( List::create(Named("V")= V_new ,
                         Named("eval_densV")=eval_densV)) ;


  }else{

    // V must NOT be sampled - log full conditional is evaluated at the current V

    // logarithm of the evaluated the density in the current point

    eval_densV = std::log(densV[ V-1 ]);

    // ENHANCED BUG DETECTION CHECK: Open segment with zero observations must have V = 0
    // Check both the extracted count and the actual count
    if(open_segment == true && (count == 0.0 || actual_count == 0) && V != 0){
      Rcpp::stop("BUG DETECTED in full_conditional_V: Open segment with zero observations has non-zero V state. V=%d, extracted_count=%.0f, actual_count=%d, start_idx=%d, end_idx=%d, LB=%.6f, UB=%.6f, end_point=%.6f, T_seg_size=%d", 
                 V, count, actual_count, start_idx_extract, end_idx_extract, LB, UB, end_point, T_seg.size());
    }

    // return a list with the current Q and the value of the full conditional in this value

    return( List::create(Named("V")=V,
                         Named("eval_densV")=eval_densV)) ;

  }


}
