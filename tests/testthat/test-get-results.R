# Tests for get_results function

context("get_results")

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

test_that("get_results returns correct structure when called with SMC_turcotte_cpp output", {
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
  expect_true("UI_bounds" %in% names(smc_result))
  expect_true("n_UI" %in% names(smc_result))
  
  # Find a non-empty Update Interval to test
  non_empty_UI <- which(!sapply(smc_result$storage_B, is.null))
  
  if (length(non_empty_UI) == 0) {
    skip("No non-empty Update Intervals found in SMC output")
  }
  
  # Test get_results on the first non-empty UI
  ui_index <- non_empty_UI[1]
  start_point <- smc_result$UI_bounds[ui_index]
  end_point <- smc_result$UI_bounds[ui_index + 1]
  
  # Call get_results
  result <- get_results(
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
  
  # Check that result is a list
  expect_type(result, "list")
  
  # Check that all required elements exist
  required_elements <- c(
    "state_container_V",
    "state_container_Z",
    "state_container_Q",
    "state_container_F",
    "num_discr_intervals",
    "discr_points"
  )
  
  for (elem in required_elements) {
    expect_true(elem %in% names(result), 
                info = paste("Missing element:", elem))
  }
  
  # Check types of returned elements
  expect_true(is.matrix(result$state_container_V))
  expect_true(is.matrix(result$state_container_Z))
  expect_true(is.matrix(result$state_container_Q))
  expect_true(is.matrix(result$state_container_F))
  expect_true(is.numeric(result$num_discr_intervals))
  expect_true(is.numeric(result$discr_points))
  
  # Check dimensions
  expect_equal(nrow(result$state_container_V), result$num_discr_intervals)
  expect_equal(nrow(result$state_container_Z), result$num_discr_intervals)
  expect_equal(nrow(result$state_container_Q), result$num_discr_intervals)
  expect_equal(nrow(result$state_container_F), result$num_discr_intervals)
  
  expect_equal(ncol(result$state_container_V), n_particle)
  expect_equal(ncol(result$state_container_Z), n_particle)
  expect_equal(ncol(result$state_container_Q), n_particle)
  expect_equal(ncol(result$state_container_F), n_particle)
  
  expect_equal(length(result$discr_points), result$num_discr_intervals)
  
  # Check that discretization points are within the interval
  expect_true(all(result$discr_points >= start_point))
  expect_true(all(result$discr_points < end_point))
  
  # Check that state values are within valid ranges
  expect_true(all(result$state_container_V >= 0 & 
                  result$state_container_V <= real_data$U))
  expect_true(all(result$state_container_Z >= 0 & 
                  result$state_container_Z <= real_data$K))
  expect_true(all(result$state_container_Q >= 0 & 
                  result$state_container_Q <= real_data$W))
  expect_true(all(result$state_container_F >= 0 & 
                  result$state_container_F <= 2))
  
  # Check no NA or NaN values
  expect_false(any(is.na(result$state_container_V)))
  expect_false(any(is.na(result$state_container_Z)))
  expect_false(any(is.na(result$state_container_Q)))
  expect_false(any(is.na(result$state_container_F)))
  expect_false(any(is.na(result$discr_points)))
  expect_false(any(is.nan(result$state_container_V)))
  expect_false(any(is.nan(result$state_container_Z)))
  expect_false(any(is.nan(result$state_container_Q)))
  expect_false(any(is.nan(result$state_container_F)))
  expect_false(any(is.nan(result$discr_points)))
})

