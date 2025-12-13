#' English words dataset
#'
#' Simulation data where the segments messages represent 5 words
#' \itemize{
#' \item MOTHER
#' \item BEACH
#' \item BEAR
#' \item FATHER
#' \item SPIDER
#' }
#' @format  a list containing parameters and variables:
#' \describe{
#'   \item{Tvec}{Vector containing the time stamps}
#'   \item{Yvec}{Vector contaning the messages}
#'   \item{Bvec}{Vector with Change-points}
#'   \item{minimum_n}{minimum number of observations allowed in each non-empty segment}
#'   \item{K}{number levels of the variable Z}
#'   \item{U}{number levels of the variable V}
#'   \item{W}{number levels of the variable Q}
#'   \item{probvec_Z}{vector with probability distribution for variable Z}
#'   \item{probvec_V}{vector with probability distribution for variable V}
#'   \item{probvec_Q}{vector with probability distribution for variable Q}
#'   \item{probvec_F}{vector with probability distribution for variable F}
#'   \item{P0}{probability to observe an empty segment after a non-empty one}
#'   \item{lambdamat}{matrix U*num_logs, each row contains the probability distribution for the message for each level of variable V}
#'   \item{alphavec}{Shape parameters ruling the rate of occurrence for the non-empty segments, for each level of variable Q}
#'   \item{muvec}{shape parameters ruling the rate of occurrence for the non-empty segments, for each level of variable Q}
#'   \item{keyvec}{shape parameters ruling the rate of occurrence for the non-empty segments, for each level of variable Z}
#'   \item{etavec}{mean parameters ruling the rate of occurrence for the non-empty segments, for each level of variable Z}
#'   \item{key0vec}{shape parameters ruling the rate of occurrence for the non-empty segments, for each level of variable F}
#'   \item{eta0vec}{mean parameters ruling the rate of occurrence for the non-empty segments, for each level of variable F}
#'   \item{num_logs}{total number of observable messages}
#'   \item{Vvec}{vector with  true value for variable V, in each segment}
#'   \item{Zvec}{vector with  true value for variable Z, in each segment}
#'   \item{Qvec}{vector with  true value for variable Q, in each segment}
#'   \item{Fvec}{vector with  true value for variable F, in each segment}
#'   \item{Mvec}{vector with  true rate in each segment}
#'   \item{Nvec}{vector with  true message count in each segment}
#' }
"english_words"

























