#' Function for extrapolating the results
#' In this function the approximated posterior distributions are extracted.
#' The observational time interval is discretized either by regularly-spaced time
#' points (\code{interval_length}) or by an explicit vector \code{discretization_points}.
#' For each discretization time the posterior distribution of the dynamic states
#' V, Z, Q and F is computed.
#' @param start_point starting point of the update interval
#' @param end_point ending point of the update interval
#' @param num_particles the number of particles generated in SMC
#' @param interval_length spacing between discretization points when
#'   \code{discretization_points} is \code{NULL}; ignored when \code{discretization_points}
#'   is supplied (breakpoint bins use consecutive grid points and \code{end_point}).
#' @param breakpoints_list list with inferred breakpoints inside the update interval, for all particles
#' @param V_state_list list with inferred V state at each breakpoint, for all particles
#' @param Z_state_list list with inferred Z state at each breakpoint, for all particles
#' @param Q_state_list list with inferred Q state at each breakpoint, for all particles
#' @param F_state_list list with inferred F state at each breakpoint, for all particles
#' @param num_states_V number of states for V
#' @param num_states_Z number of states for Z
#' @param num_states_Q number of states for Q
#' @param num_states_F number of states for F
#' @param discretization_points Optional numeric vector of discretization times in
#'   \code{[start_point, end_point)} (strict upper bound). When non-\code{NULL},
#'   these replace the regular grid from \code{seq}; must be strictly increasing
#'   after sorting and duplicate removal.
#' @export

get_results <- function(start_point,
                        end_point,
                        num_particles,
                        interval_length,
                        breakpoints_list,
                        V_state_list,
                        Z_state_list,
                        Q_state_list,
                        F_state_list,
                        num_states_V,
                        num_states_Z,
                        num_states_Q,
                        num_states_F,
                        discretization_points = NULL) {

  use_custom <- !is.null(discretization_points)
  if (use_custom) {
    d <- as.numeric(discretization_points)
    if (length(d) < 1L) {
      stop("get_results: discretization_points is empty", call. = FALSE)
    }
    if (anyNA(d) || !all(is.finite(d))) {
      stop("get_results: discretization_points must be finite and non-NA", call. = FALSE)
    }
    discr_points <- sort(unique(d))
    if (any(discr_points < start_point | discr_points >= end_point)) {
      stop(
        "get_results: each discretization point must satisfy start_point <= t < end_point",
        call. = FALSE
      )
    }
    interval_end <- c(discr_points[-1L], end_point)
  } else {
    discr_points <- seq(start_point,
                        end_point - interval_length,
                        interval_length)
    interval_end <- discr_points + interval_length
  }

  # number of discretization points created
  num_discr_intervals <-  length(discr_points)

  # initialize container matrices for storing the inferred states at each discretization time point, for each particle

  state_container_V <- matrix(0,num_discr_intervals, num_particles)
  state_container_Z <- matrix(0,num_discr_intervals, num_particles)
  state_container_Q <- matrix(0,num_discr_intervals, num_particles)
  state_container_F <- matrix(0,num_discr_intervals, num_particles)
  state_container_B <- matrix(0,num_discr_intervals, num_particles)


  for(i in 1:num_particles){

    #vector of breakpoint for particle "i"

    if(is.unsorted(unlist(breakpoints_list[i]))){

      print(paste0("start point:" , start_point))
      print(paste0("end point:" , end_point))

      stop("get_result:breakpoint vector unsorted")
    }
    # vector with segment indexes each discretization point has fallen in
    discr_point_segments <- findInterval(vec = unlist(breakpoints_list[i]),
                                   x = discr_points
                                   )


    # state dynamics at each discretization time points for particle "i"
    state_container_V[,i] <- unlist(V_state_list[i])[discr_point_segments]
    state_container_Z[,i] <- unlist(Z_state_list[i])[discr_point_segments]
    state_container_Q[,i] <- unlist(Q_state_list[i])[discr_point_segments]
    state_container_F[,i] <- unlist(F_state_list[i])[discr_point_segments]

    # Breakpoint flag per discretization bin: [discr_points[j], interval_end[j])
    # Regular grid: interval_end[j] = discr_points[j] + interval_length
    # Custom grid: interval_end[j] = discr_points[j+1] or end_point for the last bin
    breakpoints <- unlist(breakpoints_list[i])
    for (j in seq_len(num_discr_intervals)) {
      if (any(breakpoints >= discr_points[j] & breakpoints < interval_end[j])) {
        state_container_B[j, i] <- 1
      }
    }

  }

  return(list(state_container_V = state_container_V,
              state_container_Z = state_container_Z,
              state_container_Q = state_container_Q,
              state_container_F = state_container_F,
              state_container_B = state_container_B,
              num_discr_intervals = num_discr_intervals,
              discr_points = discr_points))


}
