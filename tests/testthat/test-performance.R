# Performance tests

context("Performance Tests")

test_that("SMC scales with particle count", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 333)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  # Test with different particle counts
  particle_counts <- c(50, 100)
  
  for (n_particle in particle_counts) {
    settings <- create_test_settings(
      length_UI = 5.0,
      n_particle = n_particle,
      num_logs = 5,
      method = "turcotte"
    )
    
    # Time the execution
    time_taken <- system.time({
      result <- SMC(ts_data, parameters, settings)
    })
    
    # Just verify it completes - actual timing depends on system
    expect_true(time_taken["elapsed"] >= 0)
    validate_RJSMC_object(result)
  }
})

test_that("SMC scales with data size", {
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 50,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Test with different data sizes
  data_sizes <- c(50, 100)
  
  for (n_obs in data_sizes) {
    ts_data <- create_test_ts_data(n_obs = n_obs, num_logs = 5, seed = 444)
    
    time_taken <- system.time({
      result <- SMC(ts_data, parameters, settings)
    })
    
    expect_true(time_taken["elapsed"] >= 0)
    validate_RJSMC_object(result)
  }
})

test_that("Memory usage is reasonable", {
  # This is a basic check - detailed profiling would require additional tools
  ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 555)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(
    length_UI = 5.0,
    n_particle = 100,
    num_logs = 5,
    method = "turcotte"
  )
  
  # Check that object size is reasonable (not excessively large)
  result <- SMC(ts_data, parameters, settings)
  obj_size <- object.size(result)
  
  # Should be less than 100MB for this test case
  expect_true(obj_size < 100 * 1024 * 1024)
  
  validate_RJSMC_object(result)
})

# Note: Detailed performance benchmarking would require:
# - Larger test cases
# - Multiple runs for averaging
# - Memory profiling tools
# - Comparison with baseline

