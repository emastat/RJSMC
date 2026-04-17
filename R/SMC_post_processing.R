#' Post-process SMC output into an RJSMC object
#'
#' Extracts marginal posteriors at a chosen time discretization from the list
#' returned by \code{\link{SMC}} without re-running the particle filter. Supply
#' either \code{interval_length} (regular grid) or \code{discretization_points}
#' (custom grid), not both.
#'
#' @param out_SMC_cpp List returned by \code{\link{SMC}}.
#' @param parameters Model parameter list, same as passed to \code{\link{SMC}}
#'   (uses \code{U}, \code{W}, \code{K}).
#' @param settings Settings list, same as passed to \code{\link{SMC}} (uses
#'   \code{n_particle}; optional \code{dir} for debug \code{saveRDS} dumps).
#' @param interval_length Spacing for a regular grid: pass a positive number and
#'   leave \code{discretization_points} as \code{NULL}. Mutually exclusive with
#'   \code{discretization_points}.
#' @param discretization_points Numeric vector of discretization times (same scale
#'   as \code{UI_bounds}); leave \code{interval_length} as \code{NULL}. Per update
#'   interval \code{[SP, EP)}, values with \code{SP <= t < EP} are passed to
#'   \code{\link{get_results}}. Mutually exclusive with \code{interval_length}.
#' @param elapsed_time Optional numeric seconds for the SMC C++ run; if
#'   \code{NULL}, slot \code{elapsed_time} is set to \code{NA_real_}.
#' @return Object of class \code{\linkS4class{RJSMC}}.
#' @seealso \code{\link{SMC}}
#' @export
smc_post_processing <- function(out_SMC_cpp,
                                parameters,
                                settings,
                                interval_length = NULL,
                                discretization_points = NULL,
                                elapsed_time = NULL) {

  use_interval <- !is.null(interval_length)
  use_discr <- !is.null(discretization_points)
  if (use_interval && use_discr) {
    stop(
      "smc_post_processing: pass exactly one of interval_length or discretization_points (not both)",
      call. = FALSE
    )
  }
  if (!use_interval && !use_discr) {
    stop(
      "smc_post_processing: pass exactly one of interval_length or discretization_points",
      call. = FALSE
    )
  }

  if (is.null(elapsed_time)) {
    elapsed_time <- NA_real_
  }

  if (use_discr) {
    dp_global <- as.numeric(discretization_points)
    if (anyNA(dp_global) || !all(is.finite(dp_global))) {
      stop("smc_post_processing: discretization_points must be finite and non-NA",
           call. = FALSE)
    }
  } else {
    if (!isTRUE(is.numeric(interval_length)) || length(interval_length) != 1L ||
        !is.finite(interval_length) || interval_length <= 0) {
      stop("smc_post_processing: interval_length must be a single finite positive number",
           call. = FALSE)
    }
  }

  n_particle <- settings$n_particle
  U <- parameters$U
  W <- parameters$W
  K <- parameters$K

  n_UI <- out_SMC_cpp$n_UI
  UI_bounds <- out_SMC_cpp$UI_bounds


  posteriors_container_V <- NULL
  posteriors_container_Z <- NULL
  posteriors_container_Q <- NULL
  posteriors_container_F <- NULL
  posteriors_container_B <- NULL
  points_container <- NULL
  UI_index_vector <- NULL

  non_empty_UI <- setdiff(
    1:(length(UI_bounds)-1),
    which(sapply(out_SMC_cpp$storage_B, is.null
                  )
          )
    )

  non_empty_UI_bounds = list()
  SP = UI_bounds[1]
  for(i in non_empty_UI){
    EP = UI_bounds[i+1]
    non_empty_UI_bounds = c(non_empty_UI_bounds,list(c(SP,EP)))
    SP = EP
  }

  for(i in 1:length(non_empty_UI)){

    non_empty_idx = non_empty_UI[[i]]

    sum_na = sum(is.na(unlist(out_SMC_cpp$storage_weight[[non_empty_idx]])))

    if(sum_na>0){stop("smc_post_processing: NA in storage_weight")}

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

    SP_i <- non_empty_UI_bounds[[i]][1]
    EP_i <- non_empty_UI_bounds[[i]][2]
    discr_i <- if (!use_discr) {
      NULL
    } else {
      d <- sort(unique(dp_global[dp_global >= SP_i & dp_global < EP_i]))
      if (length(d) < 1L) {
        stop(
          "smc_post_processing: no discretization_points in [SP, EP) for this update interval",
          call. = FALSE
        )
      }
      d
    }

    results_UI <- get_results(SP_i,
                              EP_i,
                              n_particle,
                              if (use_discr) NA_real_ else interval_length,
                              out_SMC_cpp$storage_B[[non_empty_idx]],
                              out_SMC_cpp$storage_V[[non_empty_idx]],
                              out_SMC_cpp$storage_Z[[non_empty_idx]],
                              out_SMC_cpp$storage_Q[[non_empty_idx]],
                              out_SMC_cpp$storage_F[[non_empty_idx]],
                              U,
                              W,
                              K,
                              2,
                              discretization_points = discr_i)


    temp_container_V <- compute_posterior(results_UI$num_discr_intervals,
                                          n_particle,
                                          results_UI$state_container_V,
                                          U+1,
                                          unlist(out_SMC_cpp$storage_weight[[non_empty_idx]]))

    temp_container_Z <- compute_posterior(results_UI$num_discr_intervals,
                                          n_particle,
                                          results_UI$state_container_Z,
                                          K+1,
                                          unlist(out_SMC_cpp$storage_weight[[non_empty_idx]]))

    temp_container_Q <- compute_posterior(results_UI$num_discr_intervals,
                                            n_particle,
                                            results_UI$state_container_Q,
                                            W+1,
                                            unlist(out_SMC_cpp$storage_weight[[non_empty_idx]]))

    temp_container_F <- compute_posterior(results_UI$num_discr_intervals,
                                            n_particle,
                                            results_UI$state_container_F,
                                            2+1,
                                            unlist(out_SMC_cpp$storage_weight[[non_empty_idx]]))

    temp_container_B <- compute_posterior_breakpoint(results_UI$num_discr_intervals,
                                                    n_particle,
                                                    results_UI$state_container_B,
                                                    unlist(out_SMC_cpp$storage_weight[[non_empty_idx]]))





    posteriors_container_V <- rbind(posteriors_container_V,temp_container_V)

    posteriors_container_Z <- rbind(posteriors_container_Z,temp_container_Z)

    posteriors_container_Q <- rbind(posteriors_container_Q,temp_container_Q)

    posteriors_container_F <- rbind(posteriors_container_F,temp_container_F)

    posteriors_container_B <- c(posteriors_container_B,temp_container_B)


    points_container <- c(points_container, results_UI$discr_points)

    UI_index_vector <- c(UI_index_vector, rep(non_empty_idx, length(results_UI$discr_points)))

  }


return(new("RJSMC",n_UI = n_UI,
            points_container = points_container,
            posteriors_container_V = posteriors_container_V,
            posteriors_container_Z = posteriors_container_Z,
            posteriors_container_Q = posteriors_container_Q,
            posteriors_container_F = posteriors_container_F,
            posteriors_container_B = posteriors_container_B,
            storage_B = out_SMC_cpp$storage_B,
            storage_V = out_SMC_cpp$storage_V,
            storage_weight = out_SMC_cpp$storage_weight,
            UI_index_vector = UI_index_vector,
            UI_bounds = UI_bounds,
            elapsed_time = elapsed_time
            ))

}
