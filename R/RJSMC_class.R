#' An S4 class to represent the output of the SMC function.
#'
#' @slot n_UI number of Update Intervals within which the SMC has been performed
#' @slot points_container list containing \code{n_UI} vector: each vector represent the time point where the marginal posterior distributions for the state variables have been approximated, for each different Update Interval
#' @slot posteriors_container_V  list containing \code{n_UI} matrices: each matrix represent the marginal posterior distribution of the state variable \code{V}, evaluated at each discretization point (stored in \code{point_containers}), for each different Update interval
#' @slot posteriors_container_Z  list containing \code{n_UI} matrices: each matrix represent the marginal posterior distribution of the state variable \code{Z}, evaluated at each discretization point (stored in \code{point_containers}), for each different Update interval
#' @slot posteriors_container_Q list containing \code{n_UI} matrices: each matrix represent the marginal posterior distribution of the state variable \code{Q}, evaluated at each discretization point (stored in \code{point_containers}), for each different Update interval
#' @slot posteriors_container_F - list containing \code{n_UI} matrices: each matrix represent the marginal posterior distribution of the state variable \code{F}, evaluated at each discretization point (stored in \code{point_containers}), for each different Update interval A length-one numeric vector
#' @slot posteriors_container_B numeric vector containing the marginal posterior probability that a breakpoint occurs at each discretization point (stored in \code{points_container}), concatenated across all Update Intervals
#' @slot UI_index_vector Index vector describing which interval each time points is evaluated within.
RJSMC <- setClass("RJSMC",
                    slots = list(n_UI = "integer",
                                 points_container="numeric",
                                 posteriors_container_V="matrix",
                                 posteriors_container_Z="matrix",
                                 posteriors_container_Q="matrix",
                                 posteriors_container_F="matrix",
                                 posteriors_container_B="numeric",
                                 UI_index_vector="integer")
)
