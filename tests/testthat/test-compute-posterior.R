# Tests for compute_posterior function

context("compute_posterior")

# Load real_data if available
load_real_data <- function() {
  if (exists("real_data", envir = .GlobalEnv)) {
    return(get("real_data", envir = .GlobalEnv))
  }
  
  # Try to load from package data
  tryCatch({
    data("real_data", package = "RJSMC", envir = environment())
    return(real_data)
  }, error = function(e) {
    return(NULL)
  })
}

test_that("compute_posterior returns correct structure when called with get_results output", {
  # Load real_data
  real_data <- load_real_data()
  
  if (is.null(real_data)) {
    skip("real_data not available")
  }
  
  # First, run SMC_turcotte_cpp to get actual output
  n_particle <- 1000
  
  smc_result <- SMC_turcotte_cpp(
    real_data$Yvec[1:800],
    real_data$Tvec[1:800],
    length_UI = 5,
    n_particle = n_particle,
    U = real_data$U,
    W = real_data$W,
    K = real_data$K,
    num_logs = real_data$num_logs,
    lambdamat = real_data$lambdamat,
    keyvec = real_data$keyvec,
    etavec = real_data$etavec,
    key0vec = real_data$key0vec,
    eta0vec = real_data$eta0vec,
    alphavec = real_data$alphavec,
    muvec = real_data$muvec,
    probvec_V = real_data$probvec_V,
    probvec_Z = real_data$probvec_Z,
    probvec_Q = real_data$probvec_Q,
    probvec_F = real_data$probvec_F,
    P0 = real_data$P0,
    minimum_n = real_data$minimum_n,
    Jss1 = 1/3,
    Js1s = 1/3,
    Smax = 150,
    n_ite = 15000,
    burn_in = 3000,
    thinning = 5
  )
  
  # Validate SMC output structure
  expect_type(smc_result, "list")
  expect_true("storage_B" %in% names(smc_result))
  expect_true("storage_V" %in% names(smc_result))
  expect_true("storage_Z" %in% names(smc_result))
  expect_true("storage_Q" %in% names(smc_result))
  expect_true("storage_F" %in% names(smc_result))
  expect_true("storage_weight" %in% names(smc_result))
  expect_true("UI_bounds" %in% names(smc_result))
  
  # Find a non-empty Update Interval to test
  non_empty_UI <- which(!sapply(smc_result$storage_B, is.null))
  
  if (length(non_empty_UI) == 0) {
    skip("No non-empty Update Intervals found in SMC output")
  }
  
  # Test compute_posterior on the first non-empty UI
  ui_index <- non_empty_UI[1]
  start_point <- smc_result$UI_bounds[ui_index]
  end_point <- smc_result$UI_bounds[ui_index + 1]
  
  # Get weights for this UI
  weight_vec <- unlist(smc_result$storage_weight[[ui_index]])
  
  # Validate weights
  expect_equal(length(weight_vec), n_particle)
  expect_false(any(is.na(weight_vec)))
  expect_false(any(is.nan(weight_vec)))
  expect_true(all(weight_vec >= 0))
  
  # Run get_results to get state containers (as done in SMC.R)
  results_UI <- get_results(
    start_point = start_point,
    end_point = end_point,
    num_particles = n_particle,
    interval_length = 0.01,
    breakpoints_list = smc_result$storage_B[[ui_index]],
    V_state_list = smc_result$storage_V[[ui_index]],
    Z_state_list = smc_result$storage_Z[[ui_index]],
    Q_state_list = smc_result$storage_Q[[ui_index]],
    F_state_list = smc_result$storage_F[[ui_index]],
    num_states_V = real_data$U,
    num_states_Z = real_data$K,
    num_states_Q = real_data$W,
    num_states_F = 2
  )
  
  # Validate get_results output
  expect_type(results_UI, "list")
  expect_true("state_container_V" %in% names(results_UI))
  expect_true("state_container_Z" %in% names(results_UI))
  expect_true("state_container_Q" %in% names(results_UI))
  expect_true("state_container_F" %in% names(results_UI))
  expect_true("num_discr_intervals" %in% names(results_UI))
  
  # Test compute_posterior for V states
  posterior_V <- compute_posterior(
    num_discr_intervals = results_UI$num_discr_intervals,
    num_particles = n_particle,
    state_container = results_UI$state_container_V,
    num_states = real_data$U + 1,
    weight_vec = weight_vec
  )
  
  # Check structure
  expect_true(is.matrix(posterior_V))
  expect_equal(nrow(posterior_V), results_UI$num_discr_intervals)
  expect_equal(ncol(posterior_V), real_data$U + 1)
  
  # Check for missing values
  expect_false(any(is.na(posterior_V)))
  expect_false(any(is.nan(posterior_V)))
  
  # Check that all values are non-negative (weights are accumulated)
  expect_true(all(posterior_V >= 0))
  
  # Check that state indices in state_container are valid
  expect_true(all(results_UI$state_container_V >= 0 & 
                  results_UI$state_container_V <= real_data$U))
  
  # Test compute_posterior for Z states
  posterior_Z <- compute_posterior(
    num_discr_intervals = results_UI$num_discr_intervals,
    num_particles = n_particle,
    state_container = results_UI$state_container_Z,
    num_states = real_data$K + 1,
    weight_vec = weight_vec
  )
  
  expect_true(is.matrix(posterior_Z))
  expect_equal(nrow(posterior_Z), results_UI$num_discr_intervals)
  expect_equal(ncol(posterior_Z), real_data$K + 1)
  expect_false(any(is.na(posterior_Z)))
  expect_false(any(is.nan(posterior_Z)))
  expect_true(all(posterior_Z >= 0))
  expect_true(all(results_UI$state_container_Z >= 0 & 
                  results_UI$state_container_Z <= real_data$K))
  
  # Test compute_posterior for Q states
  posterior_Q <- compute_posterior(
    num_discr_intervals = results_UI$num_discr_intervals,
    num_particles = n_particle,
    state_container = results_UI$state_container_Q,
    num_states = real_data$W + 1,
    weight_vec = weight_vec
  )
  
  expect_true(is.matrix(posterior_Q))
  expect_equal(nrow(posterior_Q), results_UI$num_discr_intervals)
  expect_equal(ncol(posterior_Q), real_data$W + 1)
  expect_false(any(is.na(posterior_Q)))
  expect_false(any(is.nan(posterior_Q)))
  expect_true(all(posterior_Q >= 0))
  expect_true(all(results_UI$state_container_Q >= 0 & 
                  results_UI$state_container_Q <= real_data$W))
  
  # Test compute_posterior for F states
  posterior_F <- compute_posterior(
    num_discr_intervals = results_UI$num_discr_intervals,
    num_particles = n_particle,
    state_container = results_UI$state_container_F,
    num_states = 2 + 1,  # F has 3 states (0, 1, 2)
    weight_vec = weight_vec
  )
  
  expect_true(is.matrix(posterior_F))
  expect_equal(nrow(posterior_F), results_UI$num_discr_intervals)
  expect_equal(ncol(posterior_F), 3)
  expect_false(any(is.na(posterior_F)))
  expect_false(any(is.nan(posterior_F)))
  expect_true(all(posterior_F >= 0))
  expect_true(all(results_UI$state_container_F >= 0 & 
                  results_UI$state_container_F <= 2))
  
  # Verify that compute_posterior correctly accumulates weights
  # For each discretization point, the sum across all states should equal
  # the sum of all particle weights
  total_weight <- sum(weight_vec)
  
  for (i in 1:results_UI$num_discr_intervals) {
    # Sum of posterior probabilities for V should equal total weight
    expect_equal(sum(posterior_V[i, ]), total_weight, 
                 tolerance = 1e-10,
                 info = paste("Row", i, "of posterior_V"))
    
    # Sum of posterior probabilities for Z should equal total weight
    expect_equal(sum(posterior_Z[i, ]), total_weight,
                 tolerance = 1e-10,
                 info = paste("Row", i, "of posterior_Z"))
    
    # Sum of posterior probabilities for Q should equal total weight
    expect_equal(sum(posterior_Q[i, ]), total_weight,
                 tolerance = 1e-10,
                 info = paste("Row", i, "of posterior_Q"))
    
    # Sum of posterior probabilities for F should equal total weight
    expect_equal(sum(posterior_F[i, ]), total_weight,
                 tolerance = 1e-10,
                 info = paste("Row", i, "of posterior_F"))
  }
  
  # Verify that the accumulation logic is correct
  # For a specific discretization point, manually compute what the posterior should be
  test_point <- min(10, results_UI$num_discr_intervals)
  
  # For V states at test_point
  manual_posterior_V <- rep(0, real_data$U + 1)
  for (particle_idx in 1:n_particle) {
    state_value <- results_UI$state_container_V[test_point, particle_idx]
    manual_posterior_V[state_value + 1] <- manual_posterior_V[state_value + 1] + weight_vec[particle_idx]
  }
  
  expect_equal(posterior_V[test_point, ], manual_posterior_V,
               tolerance = 1e-10,
               info = "Manual computation of posterior_V matches compute_posterior output")
  
  # Check that at least some states have non-zero probability
  # (This ensures the function is actually working, not just returning zeros)
  expect_true(any(posterior_V > 0))
  expect_true(any(posterior_Z > 0))
  expect_true(any(posterior_Q > 0))
  expect_true(any(posterior_F > 0))
})

test_that("compute_posterior handles edge cases correctly", {
  # Load real_data
  real_data <- load_real_data()
  
  if (is.null(real_data)) {
    skip("real_data not available")
  }
  
  # Test with minimal data
  n_particle <- 10
  
  smc_result <- SMC_turcotte_cpp(
    real_data$Yvec[100:400],
    real_data$Tvec[100:400],
    length_UI = 3,
    n_particle = n_particle,
    U = real_data$U,
    W = real_data$W,
    K = real_data$K,
    num_logs = real_data$num_logs,
    lambdamat = real_data$lambdamat,
    keyvec = real_data$keyvec,
    etavec = real_data$etavec,
    key0vec = real_data$key0vec,
    eta0vec = real_data$eta0vec,
    alphavec = real_data$alphavec,
    muvec = real_data$muvec,
    probvec_V = real_data$probvec_V,
    probvec_Z = real_data$probvec_Z,
    probvec_Q = real_data$probvec_Q,
    probvec_F = real_data$probvec_F,
    P0 = real_data$P0,
    minimum_n = real_data$minimum_n,
    Jss1 = 1/3,
    Js1s = 1/3,
    Smax = 50,
    n_ite = 1000,
    burn_in = 100,
    thinning = 5
  )
  
  non_empty_UI <- which(!sapply(smc_result$storage_B, is.null))
  
  if (length(non_empty_UI) == 0) {
    skip("No non-empty Update Intervals found")
  }
  
  ui_index <- non_empty_UI[1]
  start_point <- smc_result$UI_bounds[ui_index]
  end_point <- smc_result$UI_bounds[ui_index + 1]
  weight_vec <- unlist(smc_result$storage_weight[[ui_index]])
  
  results_UI <- get_results(
    start_point = start_point,
    end_point = end_point,
    num_particles = n_particle,
    interval_length = 0.01,
    breakpoints_list = smc_result$storage_B[[ui_index]],
    V_state_list = smc_result$storage_V[[ui_index]],
    Z_state_list = smc_result$storage_Z[[ui_index]],
    Q_state_list = smc_result$storage_Q[[ui_index]],
    F_state_list = smc_result$storage_F[[ui_index]],
    num_states_V = real_data$U,
    num_states_Z = real_data$K,
    num_states_Q = real_data$W,
    num_states_F = 2
  )
  
  # Test that compute_posterior works with small number of particles
  posterior_V <- compute_posterior(
    num_discr_intervals = results_UI$num_discr_intervals,
    num_particles = n_particle,
    state_container = results_UI$state_container_V,
    num_states = real_data$U + 1,
    weight_vec = weight_vec
  )
  
  expect_true(is.matrix(posterior_V))
  expect_false(any(is.na(posterior_V)))
  expect_false(any(is.nan(posterior_V)))
  expect_true(all(posterior_V >= 0))
  
  # Verify weight conservation
  total_weight <- sum(weight_vec)
  for (i in 1:min(5, nrow(posterior_V))) {
    expect_equal(sum(posterior_V[i, ]), total_weight, tolerance = 1e-10)
  }
})

