#' Post-process SMC output into an RJSMC object
#'
#' Extracts marginal posteriors at a chosen time discretization from the list
#' returned by \code{\link{SMC}}. Use a different \code{interval_length} without
#' re-running the particle filter.
#'
#' @param out_SMC_cpp List returned by \code{\link{SMC}}.
#' @param parameters Model parameter list, same as passed to \code{\link{SMC}}
#'   (uses \code{U}, \code{W}, \code{K}).
#' @param settings Settings list, same as passed to \code{\link{SMC}} (uses
#'   \code{n_particle}; optional \code{dir} for debug \code{saveRDS} dumps).
#' @param interval_length Spacing between discretization points passed to
#'   \code{\link{get_results}}.
#' @param elapsed_time Optional numeric seconds for the SMC C++ run; if
#'   \code{NULL}, slot \code{elapsed_time} is set to \code{NA_real_}.
#' @return Object of class \code{\linkS4class{RJSMC}}.
#' @seealso \code{\link{SMC}}
#' @export
smc_post_processing <- function(out_SMC_cpp,
                                parameters,
                                settings,
                                interval_length,
                                elapsed_time = NULL) {

  if (is.null(elapsed_time)) {
    elapsed_time <- NA_real_
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

    results_UI <- get_results(non_empty_UI_bounds[[i]][1],
                              non_empty_UI_bounds[[i]][2],
                              n_particle,
                              interval_length,
                              out_SMC_cpp$storage_B[[non_empty_idx]],
                              out_SMC_cpp$storage_V[[non_empty_idx]],
                              out_SMC_cpp$storage_Z[[non_empty_idx]],
                              out_SMC_cpp$storage_Q[[non_empty_idx]],
                              out_SMC_cpp$storage_F[[non_empty_idx]],
                              U,
                              W,
                              K,
                              2)


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
