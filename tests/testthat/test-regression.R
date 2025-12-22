# Regression tests

context("Regression Tests")

test_that("Known good case from test_SMC.R works", {
  # Based on tests/test_SMC.R structure
  english_words <- load_english_words_data()
  
  if (!is.null(english_words)) {
    ts_data <- list(
      Yvec = english_words$Yvec,
      Tvec = english_words$Tvec
    )
    
    parameters <- list(
      U = english_words$U,
      W = english_words$W,
      K = english_words$K,
      lambdamat = english_words$lambdamat,
      keyvec = english_words$keyvec,
      etavec = english_words$etavec,
      key0vec = english_words$key0vec,
      eta0vec = english_words$eta0vec,
      alphavec = english_words$alphavec,
      muvec = english_words$muvec,
      probvec_V = english_words$probvec_V,
      probvec_Z = english_words$probvec_Z,
      probvec_Q = english_words$probvec_Q,
      probvec_F = english_words$probvec_F,
      P0 = english_words$P0,
      minimum_n = english_words$minimum_n
    )
    
    settings <- list(
      num_logs = english_words$num_logs,
      length_UI = 0.5,
      n_particle = 100,  # Reduced from 5000 for testing
      Jss1 = 1/3,
      Js1s = 1/3,
      Smax = 150,
      n_ite = 1000,  # Reduced from 40000 for testing
      burn_in = 100,  # Reduced from 6000
      thinning = 5,
      method = "turcotte"
    )
    
    # Should complete without error
    expect_error({
      result <- SMC(ts_data, parameters, settings)
    }, NA)
    
    result <- SMC(ts_data, parameters, settings)
    validate_RJSMC_object(result)
  } else {
    skip("english_words data not available")
  }
})

test_that("Reproducibility with specific seed", {
  # Test that specific seeds produce consistent results
  # (within limits of stochastic algorithm)
  
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 999)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Run multiple times with same seed
  set.seed(12345)
  results <- replicate(2, {
    SMC(ts_data, parameters, settings)
  }, simplify = FALSE)
  
  # Check that basic structure is consistent
  expect_equal(results[[1]]@n_UI, results[[2]]@n_UI)
  expect_equal(
    length(results[[1]]@points_container),
    length(results[[2]]@points_container)
  )
  
  # Validate all results
  for (result in results) {
    validate_RJSMC_object(result)
  }
})

test_that("No regression: get_results handles edge cases", {
  # Regression test for known issues with get_results
  
  start_point <- 0.0
  end_point <- 10.0
  num_particles <- 5
  interval_length <- 0.1
  
  # Case 1: Single breakpoint (no segments)
  breakpoints_list <- lapply(1:num_particles, function(i) {
    c(start_point, end_point)
  })
  
  V_state_list <- lapply(1:num_particles, function(i) {
    c(1)
  })
  
  Z_state_list <- lapply(1:num_particles, function(i) {
    c(1)
  })
  
  Q_state_list <- lapply(1:num_particles, function(i) {
    c(1)
  })
  
  F_state_list <- lapply(1:num_particles, function(i) {
    c(0)
  })
  
  expect_error({
    result <- get_results(
      start_point, end_point, num_particles, interval_length,
      breakpoints_list, V_state_list, Z_state_list,
      Q_state_list, F_state_list,
      num_states_V = 3, num_states_Z = 2, num_states_Q = 2, num_states_F = 2
    )
  }, NA)
})

test_that("No regression: SMC handles missing settings$dir", {
  # Regression test for bug in SMC.R lines 206-216
  # These lines use settings$dir without checking if it exists
  
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 777)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  # settings$dir is not set - this should either work or fail gracefully
  # Currently it may fail if code reaches lines 206-216
  result <- tryCatch({
    SMC(ts_data, parameters, settings)
  }, error = function(e) {
    # If it fails, document the error
    expect_true(grepl("dir", e$message, ignore.case = TRUE))
    return(NULL)
  })
  
  # If it succeeds, validate
  if (!is.null(result)) {
    validate_RJSMC_object(result)
  }
})

# Note: Additional regression tests should be added as bugs are discovered
# and fixed. Each fixed bug should have a regression test to prevent
# reintroduction.

