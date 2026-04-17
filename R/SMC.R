#' SMC: Sequential Monte Carlo (particle filter + RJMCMC)
#'
#' Given a time series of time-stamps and messages, runs the particle filter with
#' Reversible Jump Markov Chain Monte Carlo for the state dynamics in Gramuglia
#' et al. (2020): Z, Q, V, and F. The observation period is split into update
#' intervals of length \code{length_UI}.
#'
#' This function returns only the raw C++ SMC output (particle states and weights
#' per update interval). To build discretized marginal posteriors and an
#' \code{\linkS4class{RJSMC}} object, call \code{\link{smc_post_processing}} with
#' the same \code{parameters} and \code{settings} lists and your chosen
#' \code{interval_length}.
#'
#' @param ts_data list with data:
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
#' @return Named list from the C++ SMC routine: \code{storage_B}, \code{storage_V},
#'   \code{storage_Z}, \code{storage_Q}, \code{storage_F}, \code{storage_S},
#'   \code{storage_weight}, integer \code{n_UI}, numeric \code{UI_bounds}.
#' @seealso \code{\link{smc_post_processing}}
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

  out_SMC_cpp
}
