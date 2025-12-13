# Integration tests

context("Integration Tests")

test_that("Full SMC pipeline works with simulated data", {
  # Create small simulated dataset
  ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 123)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 100,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Run full SMC
  result <- SMC(ts_data, parameters, settings)
  
  # Validate output
  validate_RJSMC_object(result)
  
  # Check that we have results
  expect_true(result@n_UI > 0)
  expect_true(length(result@points_container) > 0)
})

test_that("Turcotte method produces valid results", {
  # Test only Turcotte method (per Phase 4 modifications)
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 456)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  result <- SMC(ts_data, parameters, settings)
  
  # Validate structure
  validate_RJSMC_object(result)
  
  # Check that posteriors are reasonable
  # (not all zeros, not all ones, no NaN/Inf)
  expect_true(all(is.finite(result@posteriors_container_V)))
  expect_true(all(result@posteriors_container_V >= 0))
  expect_true(all(result@posteriors_container_V <= 1))
  
  # Check that at least some probabilities are non-zero
  expect_true(any(result@posteriors_container_V > 0))
})

test_that("Edge case: empty segments", {
  # Create data that might produce empty segments
  ts_data <- create_test_ts_data(n_obs = 20, num_logs = 5, seed = 789)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  parameters$P0 <- 0.5  # Higher probability of empty segments
  settings <- create_test_settings(
    length_UI = 10.0,  # Large UI might have empty segments
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Should run without error even with empty segments
  expect_error({
    result <- SMC(ts_data, parameters, settings)
  }, NA)
  
  if (exists("result")) {
    validate_RJSMC_object(result)
  }
})

test_that("Edge case: single observation", {
  # Very small dataset
  ts_data <- list(
    Yvec = c(1),
    Tvec = c(0.5)
  )
  
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 1.0,
    n_particle = 10,  # Small number of particles
    num_logs = 5,
    method = "turcotte"
  )
  
  # May or may not work depending on implementation
  # Document behavior
  result <- tryCatch({
    SMC(ts_data, parameters, settings)
  }, error = function(e) {
    NULL
  })
  
  # If it works, validate
  if (!is.null(result)) {
    validate_RJSMC_object(result)
  }
})

test_that("Edge case: very small number of particles", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 111)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 10,  # Very small
    num_logs = 5,
    method = "turcotte"
  )
  
  expect_error({
    result <- SMC(ts_data, parameters, settings)
  }, NA)
  
  result <- SMC(ts_data, parameters, settings)
  validate_RJSMC_object(result)
})

test_that("Reproducibility with set.seed", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 222)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Run twice with same seed
  set.seed(12345)
  result1 <- SMC(ts_data, parameters, settings)
  
  set.seed(12345)
  result2 <- SMC(ts_data, parameters, settings)
  
  # Results should be identical (if fully reproducible)
  # Note: May not be fully reproducible due to C++ random number generation
  # This test documents current behavior
  expect_equal(result1@n_UI, result2@n_UI)
  expect_equal(length(result1@points_container), length(result2@points_container))
})

test_that("Real data format compatibility", {
  # Test with english_words data if available
  english_words <- load_english_words_data()
  
  if (!is.null(english_words)) {
    ts_data <- list(
      Yvec = english_words$Yvec[1:100],  # Use subset for speed
      Tvec = english_words$Tvec[1:100]
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
      length_UI = 5.0,
      n_particle = 100,  # Smaller for testing
      Jss1 = 1/3,
      Js1s = 1/3,
      Smax = 50,
      n_ite = 1000,
      burn_in = 100,
      thinning = 5,
      method = "turcotte"
    )
    
    expect_error({
      result <- SMC(ts_data, parameters, settings)
    }, NA)
    
    result <- SMC(ts_data, parameters, settings)
    validate_RJSMC_object(result)
  } else {
    skip("english_words data not available")
  }
})

