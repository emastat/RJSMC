#' SMC Function implementing particle filter for the model in  Gramuglia et al. 2020 
#'
#' Given an time series of time-stamps and messages this function
#' run a particle filter with Reversible Jump Markov Chain Montecarlo
#' steps for approximating the posterior distribution of the following
#' state dynamics variables
#'  \itemize{
#'   \item Z, ruling the length of each non-empty segment
#'   \item Q, ruling the rate of occurrence of the messages in non-empty segments
#'   \item V, ruling the segment composition (in term of messages) in non-empty segments
#'   \item F, ruling the length of empty-segments
#'   }
#' The algorithm approximate these posterior distributions by splitting the observation 
#' time interval in sub-intervals defined by the custom parameter \code{length_UI}
#'
#' @param  ts_data list with data:
#' \itemize{
#'   \item {Yvec}{IntegerVector  with the observed messages}
#'   \item {Tvec}{NumericVector with  the  observed time-stamps}
#' }
#' @param parameters list with model parameters:
#' \itemize{
#'    \item {U}{Number of levels of variable V}
#'    \item {W}{Number of levels of variable Q}
#'    \item {K}{Number of levels of variable Z}
#'    \item {lambdamat}{matrix (U * num_logs) hosting the  probability mass for the messages in each V state} 
#'    \item {keyvec}{Shape parameters ruling the length of non-empty segments in different levels of Z state} 
#'    \item {etavec}{Mean parameters ruling the length of non-empty segments in different levels of Z state}
#'    \item {key0vec}{Shape parameters ruling the length of empty segments in different levels of F state}
#'    \item {eta0vec}{Mean parameters ruling the length of empty segments in different levels of F state}
#'    \item {alphavec}{Shape parameters ruling the rate of occurrence in different levels of Z state (non-empty segments)}
#'    \item {muvec}{Mean parameters ruling the rate of occurrence in different levels of Z state}
#'    \item {probvec_V}{probability mass for the V state}
#'    \item {probvec_Z}{probability mass for the Z state}
#'    \item {probvec_Q}{probability mass for the Q state}
#'    \item {probvec_F}{probability mass for the F state}
#'    \item {P0}{(const double&) Probability to observe an empty segment after a non-empty one}
#'    \item {minimum_n}{minimum number of observations that must be observed in a non-empty segment}
#' }
#' @param settings list with simulation settings parameters:
#' \itemize{
#'    \item {length_UI}{length of the each update interval}
#'    \item {n_particle}{number of particles to generate}
#'    \item {num_logs}{Total number of different messages that can be observed}
#'    \item {Jss1}{(const double&) Probability to propose a jump forward}
#'    \item {Js1s}{(const double&) Probability to propose a jump backward}
#'    \item {Smax}{maximum number of segments allowed inside the update interval}
#'    \item {n_ite}{(const int&), number of iteration to be performed in the RJMCMC}
#'    \item {burn_in}{(const int&) Number of iterations used as burn_in period.}
#'    \item {thinning}{(const int&) One every "thinning" iterations will be stored by the RJMCMC}
#'    \item {method}{(string}) set it to "turcotte" for applying the Turcotte method or to "waste_free" for applying the Chopin method for importance density approximation)
#'    \item {recycled_particles}{(int) number of particles to sample when performing the waste_free method. Default is 2}
#' }
#' @return The output of the function is a list with the following elements
#' \itemize{
#'  \item n_UI - number of Update Intervals within which the SMC has been performed
#'  \item points_container - vector containing all the time point (i.e. discretization points) where the marginal posterior distributions for the state variables have been approximated
#'  \item posteriors_container_V - matrix (length(points_container) X U) containing the marginal probability mass for state variable V, evaluated at each discretization point
#'  \item posteriors_container_Z - matrix (length(points_container) X K) containing the marginal probability mass for state variable Z, evaluated at each discretization point
#'  \item posteriors_container_Q - matrix (length(points_container) X W) containing the marginal probability mass for state variable Q, evaluated at each discretization point
#'  \item posteriors_container_F - matrix (length(points_container) X 2) containing the marginal probability mass for state variable F, evaluated at each discretization point
#'   \item UI_index_vector - vector with the Update Interval number each discretization point in "points_container" belongs to
#'  }
#' @export

SMC = function( ts_data, 
               parameters,
               settings){

  Yvec <- ts_data$Yvec 
  Tvec <- ts_data$Tvec

  length_UI <- settings$length_UI
  n_particle <- settings$n_particle
  num_logs <- settings$num_logs
  Jss1 <- settings$Jss1
  Js1s <- settings$Js1s
  Smax <- settings$Smax
  n_ite <- settings$n_ite
  burn_in <- settings$burn_in
  thinning <- settings$thinning
  method <- settings$method
  recycled_particles <- settings$recycled_particles
  if(is.null(recycled_particles)){recycled_particles = 2}

  U <- parameters$U
  W <- parameters$W
  K <- parameters$K
  lambdamat <- parameters$lambdamat
  keyvec <- parameters$keyvec
  etavec <- parameters$etavec
  key0vec <- parameters$key0vec
  eta0vec <- parameters$eta0vec
  alphavec <- parameters$alphavec
  muvec <- parameters$muvec
  probvec_V <- parameters$probvec_V
  probvec_Z <- parameters$probvec_Z
  probvec_Q <- parameters$probvec_Q
  probvec_F <- parameters$probvec_F
  P0 <- parameters$P0
  minimum_n <- parameters$minimum_n

  # run the Sequential Monte Carlo algorithm for the  whole
  # observation period
  # (from first to last observed time-stamp of the dataset)


  if(method == "turcotte"){

    out_SMC_cpp <- SMC_turcotte_cpp(Yvec,
                                    Tvec,
                                    length_UI,
                                    n_particle,
                                    U,
                                    W,
                                    K,
                                    num_logs,
                                    lambdamat,
                                    keyvec,
                                    etavec,
                                    key0vec,
                                    eta0vec,
                                    alphavec,
                                    muvec,
                                    probvec_V,
                                    probvec_Z,
                                    probvec_Q,
                                    probvec_F,
                                    P0,
                                    minimum_n,
                                    Jss1,
                                    Js1s,
                                    Smax,
                                    n_ite,
                                    burn_in,
                                    thinning)
  }else{

    out_SMC_cpp <- SMC_waste_free_cpp(Yvec,
                                      Tvec,
                                      length_UI,
                                      n_particle,
                                      U,
                                      W,
                                      K,
                                      num_logs,
                                      lambdamat,
                                      keyvec,
                                      etavec,
                                      key0vec,
                                      eta0vec,
                                      alphavec,
                                      muvec,
                                      probvec_V,
                                      probvec_Z,
                                      probvec_Q,
                                      probvec_F,
                                      P0,
                                      minimum_n,
                                      Jss1,
                                      Js1s,
                                      Smax,
                                      n_ite,
                                      burn_in,
                                      thinning,
                                      recycled_particles)

  }


  print( "SMC completed: extracting results")

  # from the output of the SMC extract the marginal posterior
  # distributions for
  # the states V, Z, Q, F computed at some fixed discretization
  # point along the
  # observation period. As the SMC computed the approximated
  # posterior distribution
  # one update interval at a time so the results will be presented

  #number of Update intervals
  n_UI <- out_SMC_cpp$n_UI

  #vector  with bounds between the Update Intervals
  UI_bounds <- out_SMC_cpp$UI_bounds


  posteriors_container_V <- NULL
  posteriors_container_Z <- NULL
  posteriors_container_Q <- NULL
  posteriors_container_F <- NULL
  points_container <- NULL
  UI_index_vector <- NULL

  non_empty_UI <- setdiff(
    1:(length(UI_bounds)-1),
    which(sapply(out_SMC_cpp$storage_B, is.null
                  )
          )
    )

  for(i in non_empty_UI){

    sum_na = sum(is.na(unlist(out_SMC_cpp$storage_weight[[i]])))

    if(sum_na>0){stop("SMC.R-->there are NA, checkout")}

    # Save debug data only if settings$dir is provided
    # These are temporary debugging files and should be optional
    if(!is.null(settings$dir) && settings$dir != ""){
      saveRDS(UI_bounds, file =paste0(settings$dir,"UI_bounds.rds"))
      saveRDS(n_particle, file =paste0(settings$dir,"n_particle.rds"))
      saveRDS(out_SMC_cpp$storage_B, file =paste0(settings$dir,"out_SMC_cpp$storage_B.rds"))
      saveRDS(out_SMC_cpp$storage_V, file =paste0(settings$dir,"out_SMC_cpp$storage_V.rds"))
      saveRDS(out_SMC_cpp$storage_Z, file =paste0(settings$dir,"out_SMC_cpp$storage_Z.rds"))
      saveRDS(out_SMC_cpp$storage_Q, file =paste0(settings$dir,"out_SMC_cpp$storage_Q.rds"))
      saveRDS(out_SMC_cpp$storage_F, file =paste0(settings$dir,"out_SMC_cpp$storage_F.rds"))
      saveRDS(U, file =paste0(settings$dir,"U.rds"))
      saveRDS(W, file =paste0(settings$dir,"W.rds"))
      saveRDS(K, file =paste0(settings$dir,"K.rds"))
    }

    results_UI <- get_results(UI_bounds[i],
                              UI_bounds[i+1],
                              n_particle,
                              0.01,
                              out_SMC_cpp$storage_B[[i]],
                              out_SMC_cpp$storage_V[[i]],
                              out_SMC_cpp$storage_Z[[i]],
                              out_SMC_cpp$storage_Q[[i]],
                              out_SMC_cpp$storage_F[[i]],
                              U,
                              W,
                              K,
                              2)


    #stop("SMC --> check results_UI")


    temp_container_V <- compute_posterior(results_UI$num_discr_intervals,
                                          n_particle,
                                          results_UI$state_container_V,
                                          U+1,
                                          unlist(out_SMC_cpp$storage_weight[[i]]))

    temp_container_Z <- compute_posterior(results_UI$num_discr_intervals,
                                          n_particle,
                                          results_UI$state_container_Z,
                                          K+1,
                                          unlist(out_SMC_cpp$storage_weight[[i]]))

    temp_container_Q <- compute_posterior(results_UI$num_discr_intervals,
                                            n_particle,
                                            results_UI$state_container_Q,
                                            W+1,
                                            unlist(out_SMC_cpp$storage_weight[[i]]))

    temp_container_F <- compute_posterior(results_UI$num_discr_intervals,
                                            n_particle,
                                            results_UI$state_container_F,
                                            2+1,
                                            unlist(out_SMC_cpp$storage_weight[[i]]))





    posteriors_container_V <- rbind(posteriors_container_V,temp_container_V)

    posteriors_container_Z <- rbind(posteriors_container_Z,temp_container_Z)

    posteriors_container_Q <- rbind(posteriors_container_Q,temp_container_Q)

    posteriors_container_F <- rbind(posteriors_container_F,temp_container_F)


    points_container <- c(points_container, results_UI$discr_points)

    UI_index_vector <- c(UI_index_vector, rep(i, length(results_UI$discr_points)))

  }


return(new("RJSMC",n_UI = n_UI,
            points_container = points_container,
            posteriors_container_V = posteriors_container_V,
            posteriors_container_Z = posteriors_container_Z,
            posteriors_container_Q = posteriors_container_Q,
            posteriors_container_F = posteriors_container_F,
            UI_index_vector = UI_index_vector
            ))

}
