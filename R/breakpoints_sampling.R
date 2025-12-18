#' @title  Function for generating a sample of Breakpoints inside an observational interval
#' Function for generating breakpoints from a kernel density inside an interval
#'
#' @param  start_point start point for the segment (t_star)
#' @param  end_point end point of the observational interval
#' @param breakpoint_list list with breakpoints vectors sampled during the RJSMCMC
#' @param sample_size number of samples to generate
#' @return list called output with following objects
#' \describe{
#' \item{logdens}{vector with log density value for  the generate breakpoint vector, for each sample }
#' \item{ S.samp}{list of vectors with the generated breakpoints in each sample }
#' @export
breakpoints_sampling = function(
    start_point,
    end_point,
    breakpoint_list,
    sample_size
)
{

  # unlist breakpoints in sorted vector and remove breakpoint smaller than start_point ( i.e t_star)

  copy_breakpoint_list = lapply(breakpoint_list,identity)

  breakpoints_flat = do.call(c,copy_breakpoint_list)
  sorted_breakpoints = sort(breakpoints_flat)
  sorted_breakpoints = unique(sorted_breakpoints)
  jittered_breakpoints = jitter(sorted_breakpoints,amount = 0.01)
  jittered_breakpoints_inside = jittered_breakpoints
  #jittered_breakpoints_inside = jittered_breakpoints[(jittered_breakpoints>start_point) &
  #                                          (jittered_breakpoints<end_point)]

  start_point_fix = min(jittered_breakpoints)-0.01
  end_point_fix = max(jittered_breakpoints)+ 0.01
  #and remove start_point
  Sv = jittered_breakpoints_inside - start_point_fix

  if(length(Sv)==1){Sv = c(Sv,Sv)}

  if(any(Sv <=0)){

    print("this is Sv")
    print(Sv)
    print("this is start_point_fix")
    print(start_point_fix)
    print("this is end_point_fix")
    print(end_point_fix)
    print("this is jittered_breakpoints_inside")
    print(jittered_breakpoints_inside)
    print("this is jittered_breakpoints")
    print(jittered_breakpoints)
    stop("Breakpoint_sampling.R Breakpoints with value <=0 have been found. Breakpoint must be positive")
  }

  # length of the observational interval
  delta = end_point_fix - start_point_fix

  #  apply logit transformation. Allow breakpoint to be in (-inf, +inf)
  Sv.tr = log((Sv/delta)/(1-Sv/delta))

  if(any(is.nan(Sv.tr))){

    print("this is breakpoint_list")
    print(breakpoint_list)
    print("this is Sv")
    print(Sv)
    print("this is start_point")
    print(start_point_fix)
    print("this is end_point")
    print(end_point_fix)
    print("this is jittered_breakpoints_inside")
    print(jittered_breakpoints_inside)
    print("this is delta")
    print(delta)
      
    stop("breakpoint_sampling.R --> check line 45")}

  # Check if all values are identical (or very close) - Mclust will hang on this
  Sv.tr_range <- max(Sv.tr) - min(Sv.tr)
  Sv.tr_sd <- if(length(Sv.tr) > 1) sd(Sv.tr) else 0
  all_identical <- (Sv.tr_range < 1e-10 || Sv.tr_sd < 1e-10)
  
  if(all_identical){
    # Add a small jitter to make values slightly different
    # This allows Mclust to fit without hanging
    jitter_amount <- max(0.01, abs(mean(Sv.tr)) * 0.01, 0.01)
    Sv.tr <- Sv.tr + rnorm(length(Sv.tr), mean = 0, sd = jitter_amount)
  }
  
  # Normal case: call Mclust
  mclust_error = TRUE
  mclust_attempts = 0

  while(mclust_error){
    mclust_attempts <- mclust_attempts + 1

    result <- try({
      fit = Mclust(Sv.tr)
                })

    if (class(result) == "try-error"){

      print("this is breakpoint_list")
      print(breakpoint_list)
      print("this is Sv")
      print(Sv)
      print("this is start_point")
      print(start_point_fix)
      print("this is end_point")
      print(end_point_fix)
      print("this is jittered_breakpoints_inside")
      print(jittered_breakpoints_inside)
      print("this is jittered_breakpoints")
      print(jittered_breakpoints)
      print("this is delta")
      print(delta)
      print("this is Sv.tr")
      print(Sv.tr)

      stop("breakpoint_sampling.R --> check line 57")

    }else{

      # fit MVN mixture to the breakpoints
      fit = Mclust(Sv.tr, verbose = FALSE)
      mclust_error = FALSE

    }
  }

  # simulate number of breakpoints from a negative binomial, based on the empirical counts
  # of breakpoints in the particles

  if(length(copy_breakpoint_list)<2){

    empirical_breakpoint_dist = rep(length(copy_breakpoint_list[[1]]),2)

  }else{

    empirical_breakpoint_dist = sapply(breakpoint_list, length)
  }

  empirical_breakpoint_dist = as.numeric(empirical_breakpoint_dist)
  unique_counts = sort(unique(empirical_breakpoint_dist))

  # Strategy for N.samp:
  # 1) If distribution is fully degenerate (single unique count), skip fitdist()
  #    and always reuse that count.
  # 2) Otherwise, first try fitting a negative binomial. If that fails (error or
  #    non-finite estimates), or if it takes more than 15 seconds, fall back to
  #    a simple discrete uniform on the observed support of counts.
  fit_N_ok   = FALSE
  use_uniform = FALSE

  if(length(unique_counts) == 1){

    # Degenerate case: all particles have the same number of breakpoints
    k = unique_counts[1]
    N.samp = rep(k, sample_size)

  }else{

    # Attempt fitdist with 15-second timeout
    fit_try <- NULL
    fitdist_timed_out <- FALSE
    
    start_time <- Sys.time()
    
    fit_try <- tryCatch({
      # Set time limit for this expression
      setTimeLimit(elapsed = 15, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
      
      # Suppress progress bar using capture.output (per Stack Overflow solution)
      # Assign result inside capture.output, then wrap in invisible to discard output
      fit_result <- NULL
      invisible(capture.output(
        fit_result <- suppressWarnings(fitdist(empirical_breakpoint_dist, "nbinom"))
      ))
      fit_result
    }, error = function(e) {
      elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
      
      # Check if this was a timeout
      if(elapsed >= 15 || grepl("time", tolower(conditionMessage(e)), fixed = TRUE)){
        fitdist_timed_out <<- TRUE
        
        # Print diagnostic information for timeout
        n_obs <- length(empirical_breakpoint_dist)
        prop_unique <- length(unique_counts) / n_obs
        count_table <- table(empirical_breakpoint_dist)
        min_count_freq <- min(count_table)
        
        cat("\n=== fitdist TIMEOUT (>15 seconds) ===\n")
        cat(sprintf("start_point: %.6f\n", start_point))
        cat(sprintf("end_point: %.6f\n", end_point))
        cat(sprintf("interval_length: %.6f\n", end_point - start_point))
        cat(sprintf("n_particles: %d\n", length(breakpoint_list)))
        cat(sprintf("empirical_breakpoint_dist length: %d\n", n_obs))
        cat(sprintf("unique_counts: %s\n", paste(unique_counts, collapse = ", ")))
        cat(sprintf("n_unique: %d\n", length(unique_counts)))
        cat(sprintf("proportion_unique: %.4f\n", prop_unique))
        cat(sprintf("count_frequencies: %s\n", paste(count_table, collapse = ", ")))
        cat(sprintf("min_frequency: %d\n", min_count_freq))
        cat(sprintf("elapsed_time: %.2f seconds\n", elapsed))
        cat("Falling back to uniform sampling on unique counts.\n")
        cat("==========================================\n\n")
      }
      
      return(structure(list(), class = "try-error"))
    }, finally = {
      # Reset time limit
      try(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), silent = TRUE)
    })

    # Check if fitdist succeeded
    fitdist_succeeded <- !fitdist_timed_out && 
                         !inherits(fit_try, "try-error") && 
                         !is.null(fit_try) &&
                         !is.null(fit_try$estimate) &&
                         all(is.finite(fit_try$estimate))
    
    if(!fitdist_succeeded){

      # fitdist timed out, failed, or produced non-finite estimates: fall back to uniform
      use_uniform = TRUE
      support <- unique_counts
      N.samp <- sample(support, size = sample_size, replace = TRUE)

    }else{

      fit.N   <- fit_try
      fit_N_ok = TRUE
      N.samp  <- rnbinom(sample_size,
                         size = fit.N$estimate[1],
                         mu   = fit.N$estimate[2])
    }
  }

  #Simulate breakpoint positions for each the N.samp breakpoints simulated, forall the "sample_size" samples & compute the overall density

  S.samp  = vector(mode='list', length=length(N.samp))

  logdens = rep(NA,sample_size)

  for(i in 1:sample_size){

    if(fit_N_ok){

      # Proper NB fit available
      logdens[i] = dnbinom(N.samp[i],
                           size=fit.N$estimate[1],
                           mu=fit.N$estimate[2],
                         log = TRUE)

    }else if(use_uniform){

      # Uniform on the observed support of counts
      logdens[i] = -log(length(unique_counts))

    }else{

      # Fully degenerate case: N.samp[i] fixed with probability 1
      logdens[i] = 0
    }

    if(N.samp[i]>0){

      #Simulate  breakpoints on transformed scale

      #foo = sort(sim(fit$modelName,fit$parameters,n=N.samp[i])[,2])

      foo <- tryCatch({
      sort(sim(fit$modelName,fit$parameters,n=N.samp[i])[,2])
      },
      error = function(e) {
        # Code to handle the error
        stop("breakpoints_sampling.R-->: ", e$message)
      })


      #Transform back
      foo_exp = delta*exp(foo)/(1+exp(foo))

      #Calculate Jacobian on log-scale
      logJac = log(foo_exp)+log(delta-foo_exp)-log(delta)

      # update the log density
      logdens[i] = logdens[i] +
        sum(dens(foo,fit$modelName,fit$parameters,log=TRUE)+logJac)

      # numerical stability
      if (logdens[i]<-100) {

        logdens[i] = -100

      }

      S.samp[[i]] = unique(c(0, foo_exp, end_point_fix-start_point_fix))

      # START CODE TO DELETE
      if (any(is.nan(S.samp[[i]]))) {
        stop("breakpoints_sampling.R --> NaN value detected in the vector.")
      }
      #END CODE TO DELETE

    }else{

      foo = c(0, end_point_fix-start_point_fix)
      S.samp[[i]] = foo



    }
  }

  # add to each sampled vector start_point

  S.samp = lapply(S.samp, function(x) x+ start_point)

  # return list with sampled breakpoints and logdensity value

  output = list()

  output$log_density = logdens
  output$S.samp = S.samp

  return(output)


}
