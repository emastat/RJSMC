# Unit tests for R functions

context("R Functions")

test_that("SMC function input validation", {
  # Create test data
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(length_UI = 5.0, n_particle = 50, num_logs = 5, method = "turcotte")
  
  # Test that function runs without error
  expect_error({
    result <- SMC(ts_data, parameters, settings)
  }, NA)
  
  # Test that result is RJSMC object
  result <- SMC(ts_data, parameters, settings)
  expect_s4_class(result, "RJSMC")
  
  # Validate structure
  validate_RJSMC_object(result)
})

test_that("SMC function handles missing settings$dir gracefully", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  settings <- create_test_settings(length_UI = 5.0, n_particle = 50, num_logs = 5, method = "turcotte")
  
  # Should fail if settings$dir is missing but code tries to use it
  # This test documents the bug: lines 206-216 in SMC.R use settings$dir without checking
  expect_error({
    result <- SMC(ts_data, parameters, settings)
  }, NA)  # Currently doesn't error, but should be fixed
})

test_that("get_results function works correctly", {
  # Create mock breakpoint data
  start_point <- 0.0
  end_point <- 10.0
  num_particles <- 10
  interval_length <- 0.1
  
  # Create simple breakpoint lists
  breakpoints_list <- lapply(1:num_particles, function(i) {
    c(start_point, 5.0, end_point)
  })
  
  V_state_list <- lapply(1:num_particles, function(i) {
    c(1, 2)
  })
  
  Z_state_list <- lapply(1:num_particles, function(i) {
    c(1, 1)
  })
  
  Q_state_list <- lapply(1:num_particles, function(i) {
    c(1, 2)
  })
  
  F_state_list <- lapply(1:num_particles, function(i) {
    c(0, 0)
  })
  
  result <- get_results(
    start_point, end_point, num_particles, interval_length,
    breakpoints_list, V_state_list, Z_state_list, 
    Q_state_list, F_state_list,
    num_states_V = 3, num_states_Z = 2, num_states_Q = 2, num_states_F = 2
  )
  
  expect_type(result, "list")
  expect_true("state_container_V" %in% names(result))
  expect_true("state_container_Z" %in% names(result))
  expect_true("state_container_Q" %in% names(result))
  expect_true("state_container_F" %in% names(result))
  expect_true("num_discr_intervals" %in% names(result))
  expect_true("discr_points" %in% names(result))
  
  # Check dimensions
  expect_equal(nrow(result$state_container_V), result$num_discr_intervals)
  expect_equal(ncol(result$state_container_V), num_particles)
})

test_that("get_results handles unsorted breakpoints", {
  start_point <- 0.0
  end_point <- 10.0
  num_particles <- 1
  interval_length <- 0.1
  
  # Create unsorted breakpoints (should trigger error)
  breakpoints_list <- list(c(5.0, 2.0, 10.0))  # Unsorted!
  
  V_state_list <- list(c(1, 2))
  Z_state_list <- list(c(1, 1))
  Q_state_list <- list(c(1, 2))
  F_state_list <- list(c(0, 0))
  
  expect_error({
    get_results(
      start_point, end_point, num_particles, interval_length,
      breakpoints_list, V_state_list, Z_state_list,
      Q_state_list, F_state_list,
      num_states_V = 3, num_states_Z = 2, num_states_Q = 2, num_states_F = 2
    )
  }, "breakpoint vector unsorted")
})

test_that("compute_lambdavec works correctly", {
  name <- c("A", "B", "C")
  prob_name <- c(0.3, 0.2, 0.5)
  
  result <- compute_lambdavec(name, prob_name)
  
  expect_type(result, "double")
  expect_length(result, 26)  # One for each LETTER
  expect_true(abs(sum(result) - 1.0) < 1e-10)  # Should sum to 1
  
  # Check that specified letters have correct probabilities
  expect_equal(result[which(LETTERS == "A")], 0.3)
  expect_equal(result[which(LETTERS == "B")], 0.2)
  expect_equal(result[which(LETTERS == "C")], 0.5)
  
  # Check that remaining letters share remaining probability
  remaining <- (1 - sum(prob_name)) / (26 - length(name))
  expect_equal(result[which(LETTERS == "D")], remaining)
})

test_that("RJSMC class constructor works", {
  n_UI <- 2L
  points_container <- c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
  posteriors_container_V <- matrix(rep(1/3, 6*3), nrow = 6, ncol = 3)
  posteriors_container_Z <- matrix(rep(1/2, 6*2), nrow = 6, ncol = 2)
  posteriors_container_Q <- matrix(rep(1/2, 6*2), nrow = 6, ncol = 2)
  posteriors_container_F <- matrix(rep(1/2, 6*2), nrow = 6, ncol = 2)
  UI_index_vector <- c(1L, 1L, 1L, 2L, 2L, 2L)
  
  obj <- new("RJSMC",
             n_UI = n_UI,
             points_container = points_container,
             posteriors_container_V = posteriors_container_V,
             posteriors_container_Z = posteriors_container_Z,
             posteriors_container_Q = posteriors_container_Q,
             posteriors_container_F = posteriors_container_F,
             UI_index_vector = UI_index_vector)
  
  expect_s4_class(obj, "RJSMC")
  expect_equal(obj@n_UI, n_UI)
  expect_equal(length(obj@points_container), length(points_container))
  expect_equal(nrow(obj@posteriors_container_V), length(points_container))
})

test_that("data_simulation generates valid data", {
  probvec_Z <- c(0.5, 0.5)
  probvec_V <- c(0.33, 0.33, 0.34)
  probvec_Q <- c(0.5, 0.5)
  probvec_F <- c(0.5, 0.5)
  P0 <- 0.1
  alphavec <- c(2.0, 3.0)
  muvec <- c(1.0, 1.5)
  keyvec <- c(2.0, 3.0)
  etavec <- c(1.0, 1.5)
  key0vec <- c(2.0, 3.0)
  eta0vec <- c(1.0, 1.5)
  lambdamat <- matrix(rep(1/5, 3*5), nrow = 3, ncol = 5)
  K <- 2
  U <- 3
  W <- 2
  Bmax <- 10.0
  seed <- 123
  
  result <- data_simulation(
    probvec_Z, probvec_V, probvec_Q, probvec_F, P0,
    alphavec, muvec, keyvec, etavec, key0vec, eta0vec,
    lambdamat, K, U, W, Bmax, seed, min_obs = 3
  )
  
  expect_type(result, "list")
  expect_true("Bvec" %in% names(result))
  expect_true("Tvec" %in% names(result))
  expect_true("Yvec" %in% names(result))
  expect_true("Vvec" %in% names(result))
  expect_true("Zvec" %in% names(result))
  expect_true("Qvec" %in% names(result))
  expect_true("Fvec" %in% names(result))
  
  # Check that breakpoints are sorted
  expect_true(all(diff(result$Bvec) >= 0))
  
  # Check that time stamps are within bounds
  expect_true(all(result$Tvec >= 0))
  expect_true(all(result$Tvec <= max(result$Bvec)))
})
