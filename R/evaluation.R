#' @title breakpoint_cx: Confusion Matrix for Breakpoint Detection
#' @description Computes and visualizes a confusion matrix comparing true breakpoint locations
#' with estimated breakpoint locations from SMC algorithm results.
#' @param obj Object of class RJSMC, the output from RJSMC::SMC (must include storage_B and storage_weight)
#' @param true_Bvec Numeric vector containing true breakpoint locations. If NULL, evaluation cannot be computed.
#' @param Tvec Numeric vector containing observation timestamps
#' @param storage_weight List containing particle weights for each Update Interval (from out_SMC_cpp$storage_weight)
#' @param UI_bounds Numeric vector containing the bounds of the Update Intervals
#' @param Vvec Numeric vector containing true V state for each segment (between consecutive breakpoints)
#' @param U Integer, total number of different V states
#' @return List containing confusion matrix and visualization
#' @export

breakpoint_cx <- function(obj,
                         true_Bvec = NULL,
                         Tvec,
                         storage_weight = NULL,
                         UI_bounds,
                         Vvec = NULL,
                         U = NULL) {
  
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
  # STEP 2b: Extract and combine V states from all Update Intervals
  # ============================================================================
  cat("Step 2b: Extracting V states from all Update Intervals...\n")
  
  # Extract storage_V from object
  storage_V <- obj@storage_V
  
  # Initialize list to store V states for each particle
  state_V_list <- vector("list", length = n_particles)
  for(p in 1:n_particles) {
    state_V_list[[p]] <- numeric(0)
  }
  
  # Loop through non-empty Update Intervals
  for(ui_idx in 1:n_UI) {
    ui_position <- non_empty_UI[ui_idx]
    ui_V_states <- storage_V[[ui_position]]
    
    # Determine if this is the last Update Interval
    is_last_UI <- (ui_idx == n_UI)
    
    # Loop through particles
    for(p in 1:n_particles) {
      # Retrieve V state vector for this Update Interval and particle
      particle_V <- unlist(ui_V_states[p])
      
      if(length(particle_V) > 0) {
        # Remove the last element UNLESS this is the last Update Interval
        if(is_last_UI) {
          # Last UI: keep all elements
          extracted_V <- particle_V
        } else {
          # Not last UI: remove the last element
          if(length(particle_V) > 1) {
            extracted_V <- particle_V[-length(particle_V)]
          } else {
            extracted_V <- numeric(0)  # Remove single element
          }
        }
        
        # Append to the corresponding particle's vector
        if(length(extracted_V) > 0) {
          state_V_list[[p]] <- c(state_V_list[[p]], extracted_V)
        }
      }
    }
  }
  
  cat(sprintf("  Extracted V states for %d particles\n", n_particles))
  
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
  
  # ============================================================================
  # STEP 5: Compute true V state and probability distribution for each observation
  # ============================================================================
  true_V_state_vec <- NULL
  V_probability_matrix <- NULL
  
  if(!is.null(Vvec) && !is.null(U)) {
    cat("Step 5: Computing true V state and probability distribution for each observation...\n")
    
    n_observations <- length(Tvec)
    
    # ========================================================================
    # Step 5a: Compute true V state for each observation
    # ========================================================================
    true_V_state_vec <- numeric(n_observations)
    
    # Sort true breakpoints for easier segment finding
    sorted_true_Bvec <- sort(true_Bvec)
    n_breakpoints <- length(sorted_true_Bvec)
    expected_Vvec_length <- n_breakpoints - 1  # Vvec contains V states for segments between consecutive breakpoints
    
    # Check that Vvec has the correct length
    if(length(Vvec) != expected_Vvec_length) {
      stop(sprintf("breakpoint_cx: Vvec length (%d) does not match expected length (%d). Vvec should contain V states for segments between consecutive breakpoints (n_breakpoints - 1 = %d).",
                   length(Vvec), expected_Vvec_length, expected_Vvec_length))
    }
    
    for(i in 1:n_observations) {
      obs_time <- Tvec[i]
      
      # Find which segment the observation falls into
      # Vvec[i] corresponds to the segment between breakpoint i and breakpoint i+1
      # If observation falls exactly on a breakpoint, it belongs to the segment after
      segment_idx <- NA
      
      # Check if observation falls before first breakpoint
      if(obs_time < sorted_true_Bvec[1]) {
        stop(sprintf("breakpoint_cx: Observation %d (time=%.6f) falls before first breakpoint (%.6f). Vvec does not contain V state for this segment.",
                     i, obs_time, sorted_true_Bvec[1]))
      }
      
      # Check if observation falls after last breakpoint
      if(obs_time >= sorted_true_Bvec[n_breakpoints]) {
        stop(sprintf("breakpoint_cx: Observation %d (time=%.6f) falls after last breakpoint (%.6f). Vvec does not contain V state for this segment.",
                     i, obs_time, sorted_true_Bvec[n_breakpoints]))
      }
      
      # Find which segment between consecutive breakpoints the observation falls into
      for(bp_idx in 1:(n_breakpoints - 1)) {
        if(obs_time >= sorted_true_Bvec[bp_idx] && obs_time < sorted_true_Bvec[bp_idx + 1]) {
          segment_idx <- bp_idx
          break
        } else if(obs_time == sorted_true_Bvec[bp_idx + 1]) {
          # Exactly on breakpoint: belongs to segment after (next segment)
          if(bp_idx + 1 < n_breakpoints) {
            segment_idx <- bp_idx + 1
          } else {
            stop(sprintf("breakpoint_cx: Observation %d (time=%.6f) falls exactly on last breakpoint (%.6f). Vvec does not contain V state for segment after last breakpoint.",
                         i, obs_time, sorted_true_Bvec[n_breakpoints]))
          }
          break
        }
      }
      
      if(is.na(segment_idx)) {
        stop(sprintf("breakpoint_cx: Could not determine segment for observation %d (time=%.6f). This should not happen.",
                     i, obs_time))
      }
      
      # Retrieve V state for this segment
      true_V_state_vec[i] <- Vvec[segment_idx]
    }
    
    # ========================================================================
    # Step 5b: Compute probability distribution for simulated V state (vectorized)
    # ========================================================================
    # Create matrix: rows = observations, columns = [timestamp, V=1, V=2, ..., V=U]
    # V states are 1-indexed: 1, 2, ..., U
    V_probability_matrix <- matrix(0, nrow = n_observations, ncol = 1 + U)
    colnames(V_probability_matrix) <- c("timestamp", paste0("V_", 1:U))
    V_probability_matrix[, 1] <- Tvec  # First column is timestamp
    
    # Vectorized computation: loop through particles, but vectorize over observations
    for(j in 1:n_particles) {
      particle_bp <- particle_breakpoints_list[[j]]
      particle_V <- state_V_list[[j]]
      particle_weights <- particle_breakpoints_weights_list[[j]]
      
      if(length(particle_bp) == 0) {
        # No breakpoints: all observations fall in the single segment
        segment_indices <- rep(1, n_observations)
      } else {
        # Use findInterval to vectorize segment finding for all observations at once
        # findInterval returns the index of the rightmost breakpoint <= obs_time
        # We need to adjust: if obs_time is in [B[i], B[i+1]), segment is i
        # findInterval with rightmost.closed=FALSE gives us the index of the left breakpoint
        segment_indices <- findInterval(Tvec, particle_bp, rightmost.closed = FALSE, left.open = FALSE)
        
        # Adjust for edge cases:
        # - If obs_time < first breakpoint, segment_indices will be 0, set to 1
        # - If obs_time >= last breakpoint, segment_indices will be length(particle_bp), keep as is
        # - For observations exactly on breakpoints, findInterval gives the index of that breakpoint
        #   We want the segment starting at that breakpoint, so we use the breakpoint index directly
        segment_indices[segment_indices == 0] <- 1
        
        # For observations at or after the last breakpoint, use the last segment
        segment_indices[segment_indices > length(particle_V)] <- length(particle_V)
      }
      
      # Validate segment indices
      invalid_segments <- which(segment_indices < 1 | segment_indices > length(particle_V))
      if(length(invalid_segments) > 0) {
        stop(sprintf("breakpoint_cx: Invalid segment indices for particle %d. Invalid indices at observations: %s",
                     j, paste(invalid_segments, collapse=", ")))
      }
      
      # Extract V states and weights for all observations at once (vectorized)
      V_state_values <- particle_V[segment_indices]
      weight_values <- particle_weights[segment_indices]
      
      # Check for V=0 (empty segments) - should not happen
      empty_segments <- which(V_state_values == 0)
      if(length(empty_segments) > 0) {
        # Get the first problematic observation to identify the segment
        first_empty <- empty_segments[1]
        seg_idx <- segment_indices[first_empty]
        
        # Determine breakpoint_start and breakpoint_end for this segment
        if(length(particle_bp) == 0) {
          breakpoint_start <- -Inf
          breakpoint_end <- Inf
        } else if(seg_idx == 1) {
          # First segment: before first breakpoint
          breakpoint_start <- -Inf
          breakpoint_end <- particle_bp[1]
        } else if(seg_idx <= length(particle_bp)) {
          # Middle segment: between breakpoints
          breakpoint_start <- particle_bp[seg_idx - 1]
          breakpoint_end <- particle_bp[seg_idx]
        } else {
          # Last segment: after last breakpoint
          breakpoint_start <- particle_bp[length(particle_bp)]
          breakpoint_end <- Inf
        }
        
        stop(sprintf("breakpoint_cx: Observations fall into empty segments (V=0) for particle %d. Observations: %s. Problematic segment: [%.6f, %.6f) (segment index: %d)",
                     j, paste(empty_segments, collapse=", "), breakpoint_start, breakpoint_end, seg_idx))
      }
      
      # Validate V state values
      invalid_V <- which(V_state_values < 1 | V_state_values > U | V_state_values != floor(V_state_values))
      if(length(invalid_V) > 0) {
        stop(sprintf("breakpoint_cx: Invalid V state values for particle %d. Invalid at observations: %s",
                     j, paste(invalid_V, collapse=", ")))
      }
      
      # Accumulate weights into matrix (vectorized)
      # For each observation, add weight to the column corresponding to its V state
      for(v_state in 1:U) {
        obs_with_v <- which(V_state_values == v_state)
        if(length(obs_with_v) > 0) {
          V_probability_matrix[obs_with_v, v_state + 1] <- V_probability_matrix[obs_with_v, v_state + 1] + weight_values[obs_with_v]
        }
      }
    }
    
    # Normalize rows that don't sum to 1 (instead of raising error)
    tolerance <- 0.01
    row_sums <- rowSums(V_probability_matrix[, 2:(U+1)])  # Sum excluding timestamp column
    rows_to_normalize <- which(abs(row_sums - 1.0) > tolerance)
    
    if(length(rows_to_normalize) > 0) {
      cat(sprintf("  Normalizing %d rows that don't sum to 1 (tolerance=%.2f)\n", length(rows_to_normalize), tolerance))
      for(i in rows_to_normalize) {
        if(row_sums[i] > 0) {
          V_probability_matrix[i, 2:(U+1)] <- V_probability_matrix[i, 2:(U+1)] / row_sums[i]
        }
      }
    }
    
    cat(sprintf("  Computed V state probabilities for %d observations\n", n_observations))
    
    # ========================================================================
    # Step 5c: Compute estimated V states and confusion matrix
    # ========================================================================
    cat("Step 5c: Computing estimated V states and confusion matrix...\n")
    
    # For each observation, find the V state with the highest probability
    # V states are in columns 2 to (U+1), so we need to subtract 1 to get the V state value
    simulated_V_state_vec <- numeric(n_observations)
    for(i in 1:n_observations) {
      # Find column with maximum probability (excluding timestamp column)
      prob_row <- V_probability_matrix[i, 2:(U+1)]
      max_col_idx <- which.max(prob_row)
      # Convert column index to V state value (V states are 1-indexed: 1, 2, ..., U)
      simulated_V_state_vec[i] <- max_col_idx
    }
    
    # Compute confusion matrix for V states (U x U matrix)
    V_confusion_matrix <- matrix(0, nrow = U, ncol = U)
    rownames(V_confusion_matrix) <- paste0("True_V_", 1:U)
    colnames(V_confusion_matrix) <- paste0("Est_V_", 1:U)
    
    for(i in 1:n_observations) {
      true_V <- true_V_state_vec[i]
      est_V <- simulated_V_state_vec[i]
      
      # Validate V state values
      if(true_V < 1 || true_V > U || true_V != floor(true_V)) {
        stop(sprintf("breakpoint_cx: Invalid true V state value %.2f for observation %d. Expected integer in [1, %d].",
                     true_V, i, U))
      }
      if(est_V < 1 || est_V > U || est_V != floor(est_V)) {
        stop(sprintf("breakpoint_cx: Invalid estimated V state value %.2f for observation %d. Expected integer in [1, %d].",
                     est_V, i, U))
      }
      
      V_confusion_matrix[true_V, est_V] <- V_confusion_matrix[true_V, est_V] + 1
    }
    
    cat(sprintf("  Computed V state confusion matrix (%dx%d)\n", U, U))
    print(V_confusion_matrix)
  } else {
    simulated_V_state_vec <- NULL
    V_confusion_matrix <- NULL
  }
  
  # ============================================================================
  # STEP 6: Compute average error measures
  # ============================================================================
  cat("Step 6: Computing average error measures...\n")
  
  # Average error for breakpoint count estimation
  # Sum of off-diagonal elements in breakpoint confusion matrix / total intervals
  total_intervals <- n_intervals
  correct_breakpoint_predictions <- sum(diag(confusion_matrix))
  incorrect_breakpoint_predictions <- sum(confusion_matrix) - correct_breakpoint_predictions
  avg_error_breakpoints <- incorrect_breakpoint_predictions / total_intervals
  
  cat(sprintf("  Breakpoint estimation: %d correct, %d incorrect out of %d intervals\n",
              correct_breakpoint_predictions, incorrect_breakpoint_predictions, total_intervals))
  cat(sprintf("  Average error for breakpoint count: %.4f\n", avg_error_breakpoints))
  
  # Average error for V state estimation
  avg_error_V_state <- NULL
  if(!is.null(V_confusion_matrix)) {
    total_observations <- n_observations
    correct_V_predictions <- sum(diag(V_confusion_matrix))
    incorrect_V_predictions <- sum(V_confusion_matrix) - correct_V_predictions
    avg_error_V_state <- incorrect_V_predictions / total_observations
    
    cat(sprintf("  V state estimation: %d correct, %d incorrect out of %d observations\n",
                correct_V_predictions, incorrect_V_predictions, total_observations))
    cat(sprintf("  Average error for V state: %.4f\n", avg_error_V_state))
  }
  
  # Return results
  result <- list(
    intervals = intervals,
    true_counts = true_breakpoint_counts,
    estimated_counts = estimated_breakpoint_count_vec,
    confusion_matrix = confusion_matrix,
    particle_breakpoints_list = particle_breakpoints_list,
    particle_breakpoints_weights_list = particle_breakpoints_weights_list,
    state_V_list = state_V_list,
    true_V_state_vec = true_V_state_vec,
    V_probability_matrix = V_probability_matrix,
    simulated_V_state_vec = simulated_V_state_vec,
    V_confusion_matrix = V_confusion_matrix,
    avg_error_breakpoints = avg_error_breakpoints,
    avg_error_V_state = avg_error_V_state
  )
  
  return(result)
}
