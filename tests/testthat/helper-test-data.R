# Helper functions for creating test data

#' Create minimal test data for SMC
#' 
#' @param n_obs Number of observations
#' @param num_logs Number of message types
#' @param seed Random seed
#' @return List with Yvec and Tvec
create_test_ts_data <- function(n_obs = 100, num_logs = 5, seed = 123) {
  set.seed(seed)
  Tvec <- cumsum(rexp(n_obs, rate = 1))
  Yvec <- sample(1:num_logs, n_obs, replace = TRUE)
  list(Yvec = Yvec, Tvec = Tvec)
}

#' Create minimal parameters list for testing
#' 
#' @param U Number of V states
#' @param W Number of Q states
#' @param K Number of Z states
#' @param num_logs Number of message types
#' @return List of parameters
create_test_parameters <- function(U = 3, W = 2, K = 2, num_logs = 5) {
  list(
    U = U,
    W = W,
    K = K,
    lambdamat = matrix(rep(1/num_logs, U * num_logs), nrow = U, ncol = num_logs),
    keyvec = rep(2.0, K),
    etavec = rep(1.0, K),
    key0vec = c(2.0, 3.0),
    eta0vec = c(1.0, 1.5),
    alphavec = rep(2.0, W),
    muvec = rep(1.0, W),
    probvec_V = rep(1/U, U),
    probvec_Z = rep(1/K, K),
    probvec_Q = rep(1/W, W),
    probvec_F = c(0.5, 0.5),
    P0 = 0.1,
    minimum_n = 0
  )
}

#' Create minimal settings list for testing
#' 
#' @param length_UI Length of update interval
#' @param n_particle Number of particles
#' @param num_logs Number of message types
#' @param method SMC method ("turcotte" or "waste_free")
#' @return List of settings
create_test_settings <- function(length_UI = 5.0, 
                                 n_particle = 100, 
                                 num_logs = 5,
                                 method = "turcotte") {
  settings <- list(
    num_logs = num_logs,
    length_UI = length_UI,
    n_particle = n_particle,
    Jss1 = 1/3,
    Js1s = 1/3,
    Smax = 50,
    n_ite = 1000,
    burn_in = 100,
    thinning = 5,
    method = method
  )
  
  if (method == "waste_free") {
    settings$recycled_particles = 2
  }
  
  settings
}

#' Validate RJSMC object structure
#' 
#' @param obj RJSMC object to validate
#' @return TRUE if valid, stops with error if invalid
validate_RJSMC_object <- function(obj) {
  # Check class
  expect_s4_class(obj, "RJSMC")
  
  # Check slots exist
  expect_true(is.integer(obj@n_UI))
  expect_true(is.numeric(obj@points_container))
  expect_true(is.matrix(obj@posteriors_container_V))
  expect_true(is.matrix(obj@posteriors_container_Z))
  expect_true(is.matrix(obj@posteriors_container_Q))
  expect_true(is.matrix(obj@posteriors_container_F))
  expect_true(is.integer(obj@UI_index_vector))
  
  # Check dimensions match
  n_points <- length(obj@points_container)
  expect_equal(nrow(obj@posteriors_container_V), n_points)
  expect_equal(nrow(obj@posteriors_container_Z), n_points)
  expect_equal(nrow(obj@posteriors_container_Q), n_points)
  expect_equal(nrow(obj@posteriors_container_F), n_points)
  expect_equal(length(obj@UI_index_vector), n_points)
  
  # Check posteriors are probability distributions (rows sum to ~1)
  expect_true(all(abs(rowSums(obj@posteriors_container_V) - 1.0) < 0.01))
  expect_true(all(abs(rowSums(obj@posteriors_container_Z) - 1.0) < 0.01))
  expect_true(all(abs(rowSums(obj@posteriors_container_Q) - 1.0) < 0.01))
  expect_true(all(abs(rowSums(obj@posteriors_container_F) - 1.0) < 0.01))
  
  # Check no NA/NaN values
  expect_false(any(is.na(obj@points_container)))
  expect_false(any(is.na(obj@posteriors_container_V)))
  expect_false(any(is.na(obj@posteriors_container_Z)))
  expect_false(any(is.na(obj@posteriors_container_Q)))
  expect_false(any(is.na(obj@posteriors_container_F)))
  
  TRUE
}

#' Load english_words test data if available
#' 
#' @return english_words data or NULL if not available
load_english_words_data <- function() {
  if (exists("english_words", envir = .GlobalEnv)) {
    return(get("english_words", envir = .GlobalEnv))
  }
  
  # Try to load from package data
  tryCatch({
    data("english_words", package = "RJSMC", envir = environment())
    return(english_words)
  }, error = function(e) {
    return(NULL)
  })
}

