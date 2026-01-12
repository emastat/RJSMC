#' @title breakpoint_cx: Confusion Matrix for Breakpoint Detection
#' @description Computes and visualizes a confusion matrix comparing true breakpoint locations
#' with estimated breakpoint locations from SMC algorithm results.
#' @param obj Object of class RJSMC, the output from RJSMC::SMC (must include storage_B and storage_weight)
#' @param true_Bvec Numeric vector containing true breakpoint locations. If NULL, evaluation cannot be computed.
#' @param Tvec Numeric vector containing observation timestamps
#' @param storage_weight List containing particle weights for each Update Interval (from out_SMC_cpp$storage_weight)
#' @param UI_bounds Numeric vector containing the bounds of the Update Intervals
#' @return List containing confusion matrix and visualization
#' @export

breakpoint_cx <- function(obj,
                         true_Bvec = NULL,
                         Tvec,
                         storage_weight = NULL,
                         UI_bounds) {
  
  # Check if true_Bvec is provided
  if(is.null(true_Bvec)) {
    stop("breakpoint_cx: true_Bvec is required for evaluation. Cannot compute confusion matrix without true breakpoint locations.")
  }


  # ============================================================================
  # STEP 1: Create intervals from unique timestamps and count true breakpoints
  # ============================================================================
  cat("Step 1: Creating intervals from unique timestamps...\n")
  
  # Get unique sorted timestamps + UI_bounds
  unique_timestamps <- sort(unique(c(UI_bounds, Tvec)))
  n_intervals <- length(unique_timestamps) - 1 
  
  if(n_intervals < 1) {
    stop("breakpoint_cx: Need at least 2 unique timestamps to create intervals")
  }
  
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
  
  # ============================================================================
  # STEP 2: Extract and combine breakpoints from all Update Intervals
  # ============================================================================
  cat("Step 2: Extracting breakpoints from all Update Intervals...\n")
  
  # Extract storage_B from object
  storage_B <- obj@storage_B
  
  # Get indices of non-empty Update Intervals
  non_empty_UI <- which(!sapply(storage_B, is.null))
  n_UI <- length(non_empty_UI)
  
  if(n_UI == 0) {
    stop("breakpoint_cx: No non-empty Update Intervals found in storage_B")
  }
  
  # Get number of particles from first non-empty UI
  n_particles <- length(storage_B[[non_empty_UI[1]]])
  
  # Initialize list to store breakpoints for each particle
  particle_breakpoints_list <- vector("list", length = n_particles)
  for(p in 1:n_particles) {
    particle_breakpoints_list[[p]] <- numeric(0)
  }

  # Initialize list to store breakpoints' weights for each particle
  particle_breakpoints_weights_list <- vector("list", length = n_particles)
  for(p in 1:n_particles) {
    particle_breakpoints_weights_list[[p]] <- numeric(0)
  }
  
  # Loop through non-empty Update Intervals
  for(ui_idx in 1:n_UI) {
    ui_position <- non_empty_UI[ui_idx]
    ui_breakpoints <- storage_B[[ui_position]]

    ui_bound_left = UI_bounds[ui_position]
    ui_bound_right = UI_bounds[ui_position+1]
    
    # Determine if this is first, last, or middle UI
    is_first_UI <- (ui_idx == 1)
    is_last_UI <- (ui_idx == n_UI)
    
    # Loop through particles
    for(p in 1:n_particles) {
      # Retrieve breakpoint vector for this Update Interval and particle
      particle_bp <- unlist(ui_breakpoints[p])
      
      if(length(particle_bp) > 0) {
        # Apply removal rules based on Update Interval position
        if(is_first_UI) {
          # First UI: remove only the last element
          if(length(particle_bp) > 1) {
            extracted_bp <- particle_bp[-length(particle_bp)]
          } else {
            extracted_bp <- particle_bp  # Keep if only one element
          }
        } else if(is_last_UI) {
          # Last UI: remove only the first element
          if(length(particle_bp) > 1) {
            extracted_bp <- particle_bp[-1]
          } else {
            extracted_bp <- particle_bp  # Keep if only one element
          }
        } else {
          # Middle UIs: remove both first and last elements
          if(length(particle_bp) > 2) {
            extracted_bp <- particle_bp[-c(1, length(particle_bp))]
          } else if(length(particle_bp) == 2) {
            extracted_bp <- numeric(0)  # Both removed
          } else {
            extracted_bp <- numeric(0)  # Single element removed
          }
        }
        
        # Append to the corresponding particle's vector
        if(length(extracted_bp) > 0) {
          particle_breakpoints_list[[p]] <- c(particle_breakpoints_list[[p]], extracted_bp)
        }

        # Append to the corresponding particle's vector of weights

        particle_bp_weights = rep(unlist(storage_weight[[ui_position]])[p], length(extracted_bp))
        if(length(extracted_bp) > 0) {
          particle_breakpoints_weights_list[[p]] <- c(particle_breakpoints_weights_list[[p]], particle_bp_weights)
        }
      }
    }
  }
  
  cat(sprintf("  Extracted breakpoints for %d particles\n", n_particles))
  
  # ============================================================================
  # STEP 3: Compute probabilities for 0, 1, and 2 breakpoints per interval
  # ============================================================================
  cat("Step 3: Computing probabilities for breakpoint counts per interval...\n")
  
  # Initialize probability columns
  intervals$P_0 <- 0.0
  intervals$P_1 <- 0.0
  intervals$P_2 <- 0.0
  
  # Loop through each interval
  for(i in 1:n_intervals) {
    interval_start <- intervals$start[i]
    interval_end <- intervals$end[i]
    
    # Loop through all particles
    for(j in 1:n_particles) {
      # Get breakpoints and weights for this particle
      particle_bp <- particle_breakpoints_list[[j]]
      particle_weights <- particle_breakpoints_weights_list[[j]]
      
      # Count how many breakpoints fall inside the interval [start, end)
      bp_in_interval <- (particle_bp >= interval_start) & (particle_bp < interval_end)
      count <- sum(bp_in_interval)
      
      if(count == 0) {
        # Continue to next particle (no weight added)
        next
      } else if(count == 1) {
        # Add the weight of that breakpoint to P_1
        bp_idx <- which(bp_in_interval)[1]
        intervals$P_1[i] <- intervals$P_1[i] + particle_weights[bp_idx]
      } else if(count == 2) {
        # Add the weight of the first breakpoint to P_2
        bp_indices <- which(bp_in_interval)
        intervals$P_2[i] <- intervals$P_2[i] + particle_weights[bp_indices[1]]
      } else {
        # Error: more than 2 breakpoints in interval
        stop(sprintf("breakpoint_cx: Particle %d has %d breakpoints in interval [%.6f, %.6f). Maximum allowed is 2.",
                     j, count, interval_start, interval_end))
      }
    }
  }
  
  # Set P_0 = 1 - P_1 - P_2 for each interval
  intervals$P_0 <- 1 - intervals$P_1 - intervals$P_2
  
  cat(sprintf("  Computed probabilities for %d intervals\n", n_intervals))
  
  # ============================================================================
  # VALIDATION: Check probability constraints (with tolerance for rounding errors)
  # ============================================================================
  tolerance <- 1e-5
  
  for(i in 1:n_intervals) {
    P_0 <- intervals$P_0[i]
    P_1 <- intervals$P_1[i]
    P_2 <- intervals$P_2[i]
    
    # Clamp probabilities to [0, 1] if within tolerance (rounding errors)
    if(P_0 < 0 && abs(P_0) <= tolerance) {
      intervals$P_0[i] <- 0.0
      P_0 <- 0.0
    }
    if(P_1 < 0 && abs(P_1) <= tolerance) {
      intervals$P_1[i] <- 0.0
      P_1 <- 0.0
    }
    if(P_2 < 0 && abs(P_2) <= tolerance) {
      intervals$P_2[i] <- 0.0
      P_2 <- 0.0
    }
    if(P_0 > 1 && (P_0 - 1) <= tolerance) {
      intervals$P_0[i] <- 1.0
      P_0 <- 1.0
    }
    if(P_1 > 1 && (P_1 - 1) <= tolerance) {
      intervals$P_1[i] <- 1.0
      P_1 <- 1.0
    }
    if(P_2 > 1 && (P_2 - 1) <= tolerance) {
      intervals$P_2[i] <- 1.0
      P_2 <- 1.0
    }
    
    sum_probs <- P_0 + P_1 + P_2
    
    # Check each probability is between 0 and 1 (after clamping)
    if(P_0 < -tolerance || P_0 > 1 + tolerance) {
      stop(sprintf("breakpoint_cx: Invalid P_0 value %.10f for interval [%.6f, %.6f). Must be between 0 and 1.",
                   P_0, intervals$start[i], intervals$end[i]))
    }
    if(P_1 < -tolerance || P_1 > 1 + tolerance) {
      stop(sprintf("breakpoint_cx: Invalid P_1 value %.10f for interval [%.6f, %.6f). Must be between 0 and 1.",
                   P_1, intervals$start[i], intervals$end[i]))
    }
    if(P_2 < -tolerance || P_2 > 1 + tolerance) {
      stop(sprintf("breakpoint_cx: Invalid P_2 value %.10f for interval [%.6f, %.6f). Must be between 0 and 1.",
                   P_2, intervals$start[i], intervals$end[i]))
    }
    
    # Check sum is not greater than 1 (with tolerance)
    if(sum_probs > 1 + tolerance) {
      stop(sprintf("breakpoint_cx: Sum of probabilities (%.10f) exceeds 1 for interval [%.6f, %.6f). P_0=%.10f, P_1=%.10f, P_2=%.10f",
                   sum_probs, intervals$start[i], intervals$end[i], P_0, P_1, P_2))
    }
  }
  
  # ============================================================================
  # STEP 4: Compute estimated breakpoint counts and confusion matrix
  # ============================================================================
  cat("Step 4: Computing estimated breakpoint counts and confusion matrix...\n")
  
  # For each interval, select which of P_0, P_1, P_2 has the greatest value
  estimated_breakpoint_count_vec <- numeric(n_intervals)
  for(i in 1:n_intervals) {
    probs <- c(intervals$P_0[i], intervals$P_1[i], intervals$P_2[i])
    max_idx <- which.max(probs)
    estimated_breakpoint_count_vec[i] <- max_idx - 1  # Convert to 0, 1, or 2
  }
  
  # Compute confusion matrix (3x3: rows = true, columns = estimated)
  confusion_matrix <- matrix(0, nrow = 3, ncol = 3)
  rownames(confusion_matrix) <- c("True_0", "True_1", "True_2")
  colnames(confusion_matrix) <- c("Est_0", "Est_1", "Est_2")
  
  for(i in 1:n_intervals) {
    true_count <- true_breakpoint_counts[i]
    est_count <- estimated_breakpoint_count_vec[i]
    
    # Ensure counts are valid (0, 1, or 2)
    if(true_count < 0 || true_count > 2) {
      stop(sprintf("breakpoint_cx: Invalid true_count %d for interval %d", true_count, i))
    }
    if(est_count < 0 || est_count > 2) {
      stop(sprintf("breakpoint_cx: Invalid estimated_count %d for interval %d", est_count, i))
    }
    
    confusion_matrix[true_count + 1, est_count + 1] <- confusion_matrix[true_count + 1, est_count + 1] + 1
  }
  
  cat(sprintf("  Computed confusion matrix:\n"))
  print(confusion_matrix)
  
  # Return results
  result <- list(
    intervals = intervals,
    true_counts = true_breakpoint_counts,
    estimated_counts = estimated_breakpoint_count_vec,
    confusion_matrix = confusion_matrix,
    particle_breakpoints_list = particle_breakpoints_list,
    particle_breakpoints_weights_list = particle_breakpoints_weights_list
  )
  
  return(result)
}
