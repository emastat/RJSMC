#' @title breakpoint_cx: Confusion Matrix for Breakpoint Detection
#' @description Computes and visualizes a confusion matrix comparing true breakpoint locations
#' with estimated breakpoint locations from SMC algorithm results.
#' @param obj Object of class RJSMC, the output from RJSMC::SMC (must include storage_B and storage_weight)
#' @param true_Bvec Numeric vector containing true breakpoint locations. If NULL, evaluation cannot be computed.
#' @param Tvec Numeric vector containing observation timestamps
#' @param storage_weight List containing particle weights for each Update Interval (from out_SMC_cpp$storage_weight)
#' @return List containing confusion matrix and visualization
#' @export

breakpoint_cx <- function(obj,
                         true_Bvec = NULL,
                         Tvec,
                         storage_weight = NULL) {
  
  # Check if true_Bvec is provided
  if(is.null(true_Bvec)) {
    stop("breakpoint_cx: true_Bvec is required for evaluation. Cannot compute confusion matrix without true breakpoint locations.")
  }
  
  # Extract storage_B from object
  storage_B <- obj@storage_B
  
  # Use storage_weight from object if not provided
  if(is.null(storage_weight)) {
    storage_weight <- obj@storage_weight
  }
  
  # Get number of particles (from first non-null storage_B element)
  non_null_idx <- which(!sapply(storage_B, is.null))[1]
  if(is.na(non_null_idx)) {
    stop("breakpoint_cx: No valid breakpoint data found in storage_B")
  }
  n_particles <- length(storage_B[[non_null_idx]])
  
  # ============================================================================
  # STEP 1: Extract complete breakpoints for each particle
  # ============================================================================
  cat("Step 1: Extracting complete breakpoints for each particle...\n")
  
  # Initialize list to store complete breakpoints for each particle
  complete_breakpoints <- vector("list", length = n_particles)
  # Initialize list to store weights for each breakpoint (mirrors complete_breakpoints)
  complete_breakpoints_weights <- vector("list", length = n_particles)
  for(p in 1:n_particles) {
    complete_breakpoints[[p]] <- numeric(0)
    complete_breakpoints_weights[[p]] <- numeric(0)
  }
  
  # Get indices of non-empty Update Intervals
  non_empty_UI <- which(!sapply(storage_B, is.null))
  n_UI <- length(non_empty_UI)
  
  if(n_UI == 0) {
    stop("breakpoint_cx: No non-empty Update Intervals found")
  }
  
  # Loop through each Update Interval
  for(i in 1:n_UI) {
    ui_idx <- non_empty_UI[i]
    ui_breakpoints <- storage_B[[ui_idx]]
    
    # Get weights for this Update Interval
    ui_weights <- storage_weight[[ui_idx]]
    if(is.null(ui_weights) || length(ui_weights) != n_particles) {
      stop(sprintf("breakpoint_cx: Invalid weights for Update Interval %d", ui_idx))
    }
    
    # Loop through each particle
    for(p in 1:n_particles) {
      particle_bp <- unlist(ui_breakpoints[p])
      particle_weight <- ui_weights[p]
      
      if(length(particle_bp) > 0 && !is.na(particle_weight) && particle_weight > 0) {
        # Apply extraction rules based on Update Interval position
        if(i == 1) {
          # First UI: keep first, discard last
          if(length(particle_bp) > 1) {
            extracted_bp <- particle_bp[-length(particle_bp)]
            extracted_weights <- rep(particle_weight, length(extracted_bp))
          } else {
            extracted_bp <- particle_bp
            extracted_weights <- particle_weight
          }
        } else if(i == n_UI) {
          # Last UI: discard first, keep last
          if(length(particle_bp) > 1) {
            extracted_bp <- particle_bp[-1]
            extracted_weights <- rep(particle_weight, length(extracted_bp))
          } else {
            extracted_bp <- particle_bp
            extracted_weights <- particle_weight
          }
        } else {
          # Middle UIs: discard both first and last
          if(length(particle_bp) > 2) {
            extracted_bp <- particle_bp[-c(1, length(particle_bp))]
            extracted_weights <- rep(particle_weight, length(extracted_bp))
          } else if(length(particle_bp) == 2) {
            extracted_bp <- numeric(0)  # Both discarded
            extracted_weights <- numeric(0)
          } else {
            extracted_bp <- numeric(0)  # Single breakpoint discarded
            extracted_weights <- numeric(0)
          }
        }
        
        # Append to complete breakpoints and weights for this particle
        complete_breakpoints[[p]] <- c(complete_breakpoints[[p]], extracted_bp)
        complete_breakpoints_weights[[p]] <- c(complete_breakpoints_weights[[p]], extracted_weights)
      }
    }
  }
  
  # Sort breakpoints for each particle (and corresponding weights)
  for(p in 1:n_particles) {
    if(length(complete_breakpoints[[p]]) > 0) {
      sort_idx <- order(complete_breakpoints[[p]])
      complete_breakpoints[[p]] <- complete_breakpoints[[p]][sort_idx]
      complete_breakpoints_weights[[p]] <- complete_breakpoints_weights[[p]][sort_idx]
    }
  }
  
  cat(sprintf("  Extracted breakpoints for %d particles\n", n_particles))
  cat(sprintf("  Number of breakpoints per particle: min=%d, max=%d, mean=%.2f\n",
              min(sapply(complete_breakpoints, length)),
              max(sapply(complete_breakpoints, length)),
              mean(sapply(complete_breakpoints, length))))
  
  # ============================================================================
  # STEP 2: Create intervals from unique timestamps and count true breakpoints
  # ============================================================================
  cat("\nStep 2: Creating intervals from unique timestamps...\n")
  
  # Get unique sorted timestamps
  unique_timestamps <- sort(unique(Tvec))
  n_intervals <- length(unique_timestamps) - 1
  
  if(n_intervals < 1) {
    stop("breakpoint_cx: Need at least 2 unique timestamps to create intervals")
  }
  
  cat(sprintf("  Found %d unique timestamps, creating %d intervals\n", 
              length(unique_timestamps), n_intervals))
  
  # Create intervals: [t_i, t_{i+1}) for i = 1, ..., n_intervals
  intervals <- data.frame(
    start = unique_timestamps[1:n_intervals],
    end = unique_timestamps[2:(n_intervals + 1)]
  )
  
  # Count true breakpoints in each interval
  true_breakpoint_counts <- numeric(n_intervals)
  for(i in 1:n_intervals) {
    count <- sum(true_Bvec >= intervals$start[i] & true_Bvec < intervals$end[i])
    true_breakpoint_counts[i] <- count
    
    # Validation: should have 0, 1, or 2 breakpoints
    if(count > 2) {
      stop(sprintf("breakpoint_cx: Interval [%.6f, %.6f) contains %d true breakpoints. Maximum allowed is 2.",
                   intervals$start[i], intervals$end[i], count))
    }
  }
  
  cat(sprintf("  True breakpoint counts: 0=%d, 1=%d, 2=%d\n",
              sum(true_breakpoint_counts == 0),
              sum(true_breakpoint_counts == 1),
              sum(true_breakpoint_counts == 2)))
  
  # ============================================================================
  # STEP 3: Compute probability matrix for breakpoint counts per interval
  # ============================================================================
  cat("\nStep 3: Computing probability matrix for breakpoint counts...\n")
  
  # Initialize probability matrix: rows = intervals, columns = [prob_0, prob_1, prob_2]
  prob_matrix <- matrix(0, nrow = n_intervals, ncol = 3)
  colnames(prob_matrix) <- c("prob_0", "prob_1", "prob_2")
  
  # Get UI bounds to determine which UI contains each interval
  UI_bounds <- obj@UI_index_vector  # Actually, we need the actual UI bounds
  # We'll need to get this from the object or compute it
  # For now, let's get it from storage_B structure
  
  # Find which Update Intervals contain each interval
  # We'll use the midpoint of each interval to determine which UI it belongs to
  # Actually, we need the actual UI bounds - let's compute them from the data
  
  # For each interval, determine which Update Intervals it spans
  # Use product of weights across UIs (representing trajectory probability)
  # or use weights from the UI that contains the interval midpoint
  
  # For each interval
  for(i in 1:n_intervals) {
    interval_start <- intervals$start[i]
    interval_end <- intervals$end[i]
    
    # Initialize weight sums for 0, 1, 2 breakpoints
    weight_sum_0 <- 0
    weight_sum_1 <- 0
    weight_sum_2 <- 0
    
    # For each particle
    for(p in 1:n_particles) {
      # Get breakpoints and their weights for this particle
      particle_bp <- complete_breakpoints[[p]]
      particle_bp_weights <- complete_breakpoints_weights[[p]]
      
      # Find which breakpoints fall in this interval
      bp_in_interval_mask <- particle_bp >= interval_start & particle_bp < interval_end
      bp_in_interval_count <- sum(bp_in_interval_mask)
      
      # Validation: should have 0, 1, or 2 breakpoints
      if(bp_in_interval_count > 2) {
        stop(sprintf("breakpoint_cx: Particle %d has %d breakpoints in interval [%.6f, %.6f). Maximum allowed is 2.",
                     p, bp_in_interval_count, interval_start, interval_end))
      }
      
      # Sum weights of breakpoints in this interval
      if(bp_in_interval_count == 0) {
        # No breakpoints in interval: use particle's "base weight" 
        # We can use the average weight across all breakpoints, or use a default
        # For particles with no breakpoints, we'll use equal weight
        # Actually, for 0 breakpoints, we need to think about what weight to use
        # Let's use the average of all breakpoint weights for this particle, or 1/n_particles if no breakpoints
        if(length(particle_bp_weights) > 0) {
          particle_weight <- mean(particle_bp_weights)
        } else {
          particle_weight <- 1 / n_particles
        }
        weight_sum_0 <- weight_sum_0 + particle_weight
      } else if(bp_in_interval_count == 1) {
        # One breakpoint in interval: use its weight
        bp_weight <- particle_bp_weights[bp_in_interval_mask]
        weight_sum_1 <- weight_sum_1 + bp_weight
      } else {
        # Two breakpoints in interval: sum their weights
        bp_weights <- particle_bp_weights[bp_in_interval_mask]
        total_weight <- sum(bp_weights)
        weight_sum_2 <- weight_sum_2 + total_weight
      }
    }
    
    # Normalize to get probabilities
    total_weight <- weight_sum_0 + weight_sum_1 + weight_sum_2
    if(total_weight > 0) {
      prob_matrix[i, 1] <- weight_sum_0 / total_weight
      prob_matrix[i, 2] <- weight_sum_1 / total_weight
      prob_matrix[i, 3] <- weight_sum_2 / total_weight
    } else {
      # If no weight, assign equal probabilities
      prob_matrix[i, ] <- c(1/3, 1/3, 1/3)
    }
  }
  
  cat("  Probability matrix computed\n")
  cat(sprintf("  Mean probabilities: prob_0=%.3f, prob_1=%.3f, prob_2=%.3f\n",
              mean(prob_matrix[, 1]), mean(prob_matrix[, 2]), mean(prob_matrix[, 3])))
  
  # ============================================================================
  # STEP 4: Estimate breakpoint counts based on maximum probability
  # ============================================================================
  cat("\nStep 4: Estimating breakpoint counts...\n")
  
  estimated_breakpoint_counts <- numeric(n_intervals)
  for(i in 1:n_intervals) {
    # Find which column has maximum probability (0, 1, or 2)
    max_idx <- which.max(prob_matrix[i, ])
    estimated_breakpoint_counts[i] <- max_idx - 1  # Convert to 0, 1, or 2
  }
  
  cat(sprintf("  Estimated breakpoint counts: 0=%d, 1=%d, 2=%d\n",
              sum(estimated_breakpoint_counts == 0),
              sum(estimated_breakpoint_counts == 1),
              sum(estimated_breakpoint_counts == 2)))
  
  # ============================================================================
  # STEP 5: Compute and visualize confusion matrix
  # ============================================================================
  cat("\nStep 5: Computing confusion matrix...\n")
  
  # Create confusion matrix: rows = true, columns = estimated
  confusion_matrix <- table(
    True = factor(true_breakpoint_counts, levels = 0:2),
    Estimated = factor(estimated_breakpoint_counts, levels = 0:2)
  )
  
  # Convert to matrix for better formatting
  confusion_matrix <- as.matrix(confusion_matrix)
  rownames(confusion_matrix) <- paste0("True_", 0:2)
  colnames(confusion_matrix) <- paste0("Est_", 0:2)
  
  cat("  Confusion matrix:\n")
  print(confusion_matrix)
  
  # Compute summary statistics
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  cat(sprintf("\n  Overall Accuracy: %.4f\n", accuracy))
  
  # Return results
  result <- list(
    confusion_matrix = confusion_matrix,
    true_counts = true_breakpoint_counts,
    estimated_counts = estimated_breakpoint_counts,
    intervals = intervals,
    prob_matrix = prob_matrix,
    accuracy = accuracy,
    complete_breakpoints = complete_breakpoints,
    complete_breakpoints_weights = complete_breakpoints_weights
  )
  
  cat("\nEvaluation complete!\n")
  
  return(result)
}

