# Comprehensive tests for SMC_turcotte_cpp function

context("SMC_turcotte_cpp")

# ============================================================================
# Helper function to validate SMC_turcotte_cpp output structure
# ============================================================================
validate_SMC_turcotte_output <- function(result, n_particle) {
  # Check result is a list
  expect_type(result, "list")
  
  # Check all required elements exist
  required_elements <- c(
    "storage_B", "storage_V", "storage_Z", "storage_Q", 
    "storage_F", "storage_S", "storage_weight", "n_UI", "UI_bounds"
  )
  for (elem in required_elements) {
    expect_true(elem %in% names(result), 
                info = paste("Missing element:", elem))
  }
  
  # Check n_UI is a positive integer
  expect_type(result$n_UI, "integer")
  expect_true(result$n_UI > 0)
  expect_true(is.finite(result$n_UI))
  
  # Check UI_bounds is a numeric vector with n_UI+1 elements (all boundary points)
  expect_type(result$UI_bounds, "double")
  expect_length(result$UI_bounds, result$n_UI + 1)
  expect_true(all(is.finite(result$UI_bounds)))
  expect_true(all(diff(result$UI_bounds) > 0))  # Should be strictly increasing
  
  # Check storage containers are lists with n_UI elements
  storage_elements <- c("storage_B", "storage_V", "storage_Z", 
                        "storage_Q", "storage_F", "storage_S", "storage_weight")
  for (elem in storage_elements) {
    expect_type(result[[elem]], "list")
    expect_length(result[[elem]], result$n_UI)
  }
  
  # Identify non-empty update intervals (NULL intervals are empty by design)
  # This matches the logic in SMC.R lines 195-200
  non_empty_UI <- setdiff(
    1:result$n_UI,
    which(sapply(result$storage_B, is.null))
  )
  
  # Check that at least some intervals are non-empty
  expect_true(length(non_empty_UI) > 0, 
              info = "At least one update interval should have observations")
  
  # Check each NON-EMPTY update interval has correct structure
  # (Empty intervals are NULL by design and are skipped in SMC.R)
  for (j in non_empty_UI) {
    # storage_B: list of breakpoint vectors (one per particle)
    expect_type(result$storage_B[[j]], "list")
    expect_length(result$storage_B[[j]], n_particle)
    
    # storage_V, Z, Q, F: list of state vectors (one per particle)
    for (state in c("storage_V", "storage_Z", "storage_Q", "storage_F")) {
      expect_type(result[[state]][[j]], "list")
      expect_length(result[[state]][[j]], n_particle)
    }
    
    # storage_S: vector with number of segments per particle
    expect_type(result$storage_S[[j]], "integer")
    expect_length(result$storage_S[[j]], n_particle)
    expect_true(all(result$storage_S[[j]] >= 0))
    expect_true(all(is.finite(result$storage_S[[j]])))
    
    # storage_weight: list of weight vectors (one per particle)
    expect_type(result$storage_weight[[j]], "list")
    expect_length(result$storage_weight[[j]], n_particle)
    
    # Validate a few particles in detail
    for (p in 1:min(3, n_particle)) {
      S_p <- result$storage_S[[j]][p]
      
      # Breakpoints: Bvec_final has length S (number of segments),
      # containing the last breakpoint from the previous UI plus S-1 internal breakpoints
      Bvec_p <- result$storage_B[[j]][[p]]
      expect_true(length(Bvec_p) >= S_p)
      if (S_p > 0) {
        expect_true(all(is.finite(Bvec_p[1:S_p])))
      }
      
      # States: should have S elements (padded to Smax)
      for (state in c("storage_V", "storage_Z", "storage_Q", "storage_F")) {
        state_vec <- result[[state]][[j]][[p]]
        expect_true(length(state_vec) >= S_p)  # May be padded to Smax
        expect_type(state_vec, "integer")
        # Check meaningful part
        if (S_p > 0) {
          expect_true(all(is.finite(state_vec[1:S_p])))
        }
      }
      
      # Weights: should be numeric
      expect_type(result$storage_weight[[j]][[p]], "double")
      expect_length(result$storage_weight[[j]][[p]], 1)
      expect_true(is.finite(result$storage_weight[[j]][[p]]))
    }
  }
}

# ============================================================================
# Tests with Mock/Simulated Data
# ============================================================================

# test_that("SMC_turcotte_cpp runs without error with mock data", {
#   ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 123)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
#   expect_error({
#     result <- SMC_turcotte_cpp(
#       ts_data$Yvec, ts_data$Tvec,
#       length_UI = 5.0,
#       n_particle = 50,
#       U = parameters$U,
#       W = parameters$W,
#       K = parameters$K,
#       num_logs = 5,
#       lambdamat = parameters$lambdamat,
#       keyvec = parameters$keyvec,
#       etavec = parameters$etavec,
#       key0vec = parameters$key0vec,
#       eta0vec = parameters$eta0vec,
#       alphavec = parameters$alphavec,
#       muvec = parameters$muvec,
#       probvec_V = parameters$probvec_V,
#       probvec_Z = parameters$probvec_Z,
#       probvec_Q = parameters$probvec_Q,
#       probvec_F = parameters$probvec_F,
#       P0 = parameters$P0,
#       minimum_n = parameters$minimum_n,
#       Jss1 = 1/3,
#       Js1s = 1/3,
#       Smax = 50,
#       n_ite = 1000,
#       burn_in = 100,
#       thinning = 5
#     )
#   }, NA)
# })

# test_that("SMC_turcotte_cpp returns correct structure with mock data", {
#   ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 456)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
#   n_particle <- 50
  
#   result <- SMC_turcotte_cpp(
#     ts_data$Yvec, ts_data$Tvec,
#     length_UI = 5.0,
#     n_particle = n_particle,
#     U = parameters$U,
#     W = parameters$W,
#     K = parameters$K,
#     num_logs = 5,
#     lambdamat = parameters$lambdamat,
#     keyvec = parameters$keyvec,
#     etavec = parameters$etavec,
#     key0vec = parameters$key0vec,
#     eta0vec = parameters$eta0vec,
#     alphavec = parameters$alphavec,
#     muvec = parameters$muvec,
#     probvec_V = parameters$probvec_V,
#     probvec_Z = parameters$probvec_Z,
#     probvec_Q = parameters$probvec_Q,
#     probvec_F = parameters$probvec_F,
#     P0 = parameters$P0,
#     minimum_n = parameters$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 50,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   validate_SMC_turcotte_output(result, n_particle)
# })

# test_that("SMC_turcotte_cpp state values are within valid ranges (mock data)", {
#   ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 789)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
#   n_particle <- 30
  
#   result <- SMC_turcotte_cpp(
#     ts_data$Yvec, ts_data$Tvec,
#     length_UI = 5.0,
#     n_particle = n_particle,
#     U = parameters$U,
#     W = parameters$W,
#     K = parameters$K,
#     num_logs = 5,
#     lambdamat = parameters$lambdamat,
#     keyvec = parameters$keyvec,
#     etavec = parameters$etavec,
#     key0vec = parameters$key0vec,
#     eta0vec = parameters$eta0vec,
#     alphavec = parameters$alphavec,
#     muvec = parameters$muvec,
#     probvec_V = parameters$probvec_V,
#     probvec_Z = parameters$probvec_Z,
#     probvec_Q = parameters$probvec_Q,
#     probvec_F = parameters$probvec_F,
#     P0 = parameters$P0,
#     minimum_n = parameters$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 50,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   # Identify non-empty update intervals (NULL intervals are empty by design)
#   non_empty_UI <- setdiff(
#     1:result$n_UI,
#     which(sapply(result$storage_B, is.null))
#   )
  
#   # Check state values across non-empty update intervals and particles
#   for (j in non_empty_UI) {
#     for (p in 1:n_particle) {
#       S_p <- result$storage_S[[j]][p]
      
#       if (S_p > 0) {
#         # V states: 0 to U (0 for empty, 1 to U for non-empty)
#         Vvec <- result$storage_V[[j]][[p]]
#         expect_true(all(Vvec >= 0))
#         expect_true(all(Vvec <= parameters$U))
        
#         # Z states: 0 to K (0 for empty, 1 to K for non-empty)
#         Zvec <- result$storage_Z[[j]][[p]]
#         expect_true(all(Zvec >= 0))
#         expect_true(all(Zvec <= parameters$K))
        
#         # Q states: 0 to W (0 for empty, 1 to W for non-empty)
#         Qvec <- result$storage_Q[[j]][[p]]
#         expect_true(all(Qvec >= 0))
#         expect_true(all(Qvec <= parameters$W))
        
#         # F states: 0 to 2 (0 for non-empty, 1 to 2 for empty)
#         Fvec <- result$storage_F[[j]][[p]]
#         expect_true(all(Fvec >= 0))
#         expect_true(all(Fvec <= 2))
        
#         # Breakpoints should be sorted
#         Bvec <- result$storage_B[[j]][[p]]
#         expect_true(all(diff(Bvec) >= 0))
#       }
#     }
#   }
# })

# test_that("SMC_turcotte_cpp handles different particle counts (mock data)", {
#   ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 111)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
#   for (n_particle in c(10, 50, 100)) {
#     result <- SMC_turcotte_cpp(
#       ts_data$Yvec, ts_data$Tvec,
#       length_UI = 5.0,
#       n_particle = n_particle,
#       U = parameters$U,
#       W = parameters$W,
#       K = parameters$K,
#       num_logs = 5,
#       lambdamat = parameters$lambdamat,
#       keyvec = parameters$keyvec,
#       etavec = parameters$etavec,
#       key0vec = parameters$key0vec,
#       eta0vec = parameters$eta0vec,
#       alphavec = parameters$alphavec,
#       muvec = parameters$muvec,
#       probvec_V = parameters$probvec_V,
#       probvec_Z = parameters$probvec_Z,
#       probvec_Q = parameters$probvec_Q,
#       probvec_F = parameters$probvec_F,
#       P0 = parameters$P0,
#       minimum_n = parameters$minimum_n,
#       Jss1 = 1/3,
#       Js1s = 1/3,
#       Smax = 50,
#       n_ite = 1000,
#       burn_in = 100,
#       thinning = 5
#     )
    
#     validate_SMC_turcotte_output(result, n_particle)
#   }
# })

# test_that("SMC_turcotte_cpp produces reproducible results with same seed (mock data)", {
#   ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 333)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
#   n_particle <- 30
  
#   # Run twice with same seed
#   set.seed(12345)
#   result1 <- SMC_turcotte_cpp(
#     ts_data$Yvec, ts_data$Tvec,
#     length_UI = 5.0,
#     n_particle = n_particle,
#     U = parameters$U,
#     W = parameters$W,
#     K = parameters$K,
#     num_logs = 5,
#     lambdamat = parameters$lambdamat,
#     keyvec = parameters$keyvec,
#     etavec = parameters$etavec,
#     key0vec = parameters$key0vec,
#     eta0vec = parameters$eta0vec,
#     alphavec = parameters$alphavec,
#     muvec = parameters$muvec,
#     probvec_V = parameters$probvec_V,
#     probvec_Z = parameters$probvec_Z,
#     probvec_Q = parameters$probvec_Q,
#     probvec_F = parameters$probvec_F,
#     P0 = parameters$P0,
#     minimum_n = parameters$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 50,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   set.seed(12345)
#   result2 <- SMC_turcotte_cpp(
#     ts_data$Yvec, ts_data$Tvec,
#     length_UI = 5.0,
#     n_particle = n_particle,
#     U = parameters$U,
#     W = parameters$W,
#     K = parameters$K,
#     num_logs = 5,
#     lambdamat = parameters$lambdamat,
#     keyvec = parameters$keyvec,
#     etavec = parameters$etavec,
#     key0vec = parameters$key0vec,
#     eta0vec = parameters$eta0vec,
#     alphavec = parameters$alphavec,
#     muvec = parameters$muvec,
#     probvec_V = parameters$probvec_V,
#     probvec_Z = parameters$probvec_Z,
#     probvec_Q = parameters$probvec_Q,
#     probvec_F = parameters$probvec_F,
#     P0 = parameters$P0,
#     minimum_n = parameters$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 50,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   # Results should be identical
#   expect_equal(result1$n_UI, result2$n_UI)
#   expect_equal(result1$UI_bounds, result2$UI_bounds)
#   expect_equal(result1$storage_S, result2$storage_S)
  
#   # Check a few particles match
#   for (j in 1:min(3, result1$n_UI)) {
#     expect_equal(result1$storage_S[[j]], result2$storage_S[[j]])
#     for (p in 1:min(3, n_particle)) {
#       expect_equal(result1$storage_B[[j]][[p]], result2$storage_B[[j]][[p]])
#       expect_equal(result1$storage_V[[j]][[p]], result2$storage_V[[j]][[p]])
#     }
#   }
# })

# test_that("SMC_turcotte_cpp handles edge case with few observations (mock data)", {
#   # Very small dataset
#   ts_data <- create_test_ts_data(n_obs = 10, num_logs = 5, seed = 444)
#   parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
#   n_particle <- 20
  
#   result <- SMC_turcotte_cpp(
#     ts_data$Yvec, ts_data$Tvec,
#     length_UI = 2.0,
#     n_particle = n_particle,
#     U = parameters$U,
#     W = parameters$W,
#     K = parameters$K,
#     num_logs = 5,
#     lambdamat = parameters$lambdamat,
#     keyvec = parameters$keyvec,
#     etavec = parameters$etavec,
#     key0vec = parameters$key0vec,
#     eta0vec = parameters$eta0vec,
#     alphavec = parameters$alphavec,
#     muvec = parameters$muvec,
#     probvec_V = parameters$probvec_V,
#     probvec_Z = parameters$probvec_Z,
#     probvec_Q = parameters$probvec_Q,
#     probvec_F = parameters$probvec_F,
#     P0 = parameters$P0,
#     minimum_n = parameters$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 50,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   validate_SMC_turcotte_output(result, n_particle)
# })

# ============================================================================
# Tests with english_words Real Data
# ============================================================================

# test_that("SMC_turcotte_cpp runs without error with english_words data", {
#   english_words <- load_english_words_data()
  
#   if (is.null(english_words)) {
#     skip("english_words data not available")
#   }
  
#   expect_error({
#     result <- SMC_turcotte_cpp(
#       english_words$Yvec,
#       english_words$Tvec,
#       length_UI = 0.5,
#       n_particle = 100,
#       U = english_words$U,
#       W = english_words$W,
#       K = english_words$K,
#       num_logs = english_words$num_logs,
#       lambdamat = english_words$lambdamat,
#       keyvec = english_words$keyvec,
#       etavec = english_words$etavec,
#       key0vec = english_words$key0vec,
#       eta0vec = english_words$eta0vec,
#       alphavec = english_words$alphavec,
#       muvec = english_words$muvec,
#       probvec_V = english_words$probvec_V,
#       probvec_Z = english_words$probvec_Z,
#       probvec_Q = english_words$probvec_Q,
#       probvec_F = english_words$probvec_F,
#       P0 = english_words$P0,
#       minimum_n = english_words$minimum_n,
#       Jss1 = 1/3,
#       Js1s = 1/3,
#       Smax = 150,
#       n_ite = 1000,
#       burn_in = 100,
#       thinning = 5
#     )
#   }, NA)
# })

# test_that("SMC_turcotte_cpp returns correct structure with english_words data", {
#   english_words <- load_english_words_data()
  
#   if (is.null(english_words)) {
#     skip("english_words data not available")
#   }
  
#   n_particle <- 2000
  
#   result <- SMC_turcotte_cpp(
#     english_words$Yvec,
#     english_words$Tvec,
#     length_UI = 1.5,
#     n_particle = n_particle,
#     U = english_words$U,
#     W = english_words$W,
#     K = english_words$K,
#     num_logs = english_words$num_logs,
#     lambdamat = english_words$lambdamat,
#     keyvec = english_words$keyvec,
#     etavec = english_words$etavec,
#     key0vec = english_words$key0vec,
#     eta0vec = english_words$eta0vec,
#     alphavec = english_words$alphavec,
#     muvec = english_words$muvec,
#     probvec_V = english_words$probvec_V,
#     probvec_Z = english_words$probvec_Z,
#     probvec_Q = english_words$probvec_Q,
#     probvec_F = english_words$probvec_F,
#     P0 = english_words$P0,
#     minimum_n = english_words$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 150,
#     n_ite = 20000,
#     burn_in = 5000,
#     thinning = 5
#   )
  
#   validate_SMC_turcotte_output(result, n_particle)
# })



# test_that("SMC_turcotte_cpp produces reproducible results with english_words data", {
#   english_words <- load_english_words_data()
  
#   if (is.null(english_words)) {
#     skip("english_words data not available")
#   }
  
#   n_particle <- 50
  
#   # Run twice with same seed
#   set.seed(54321)
#   result1 <- SMC_turcotte_cpp(
#     english_words$Yvec,
#     english_words$Tvec,
#     length_UI = 0.5,
#     n_particle = n_particle,
#     U = english_words$U,
#     W = english_words$W,
#     K = english_words$K,
#     num_logs = english_words$num_logs,
#     lambdamat = english_words$lambdamat,
#     keyvec = english_words$keyvec,
#     etavec = english_words$etavec,
#     key0vec = english_words$key0vec,
#     eta0vec = english_words$eta0vec,
#     alphavec = english_words$alphavec,
#     muvec = english_words$muvec,
#     probvec_V = english_words$probvec_V,
#     probvec_Z = english_words$probvec_Z,
#     probvec_Q = english_words$probvec_Q,
#     probvec_F = english_words$probvec_F,
#     P0 = english_words$P0,
#     minimum_n = english_words$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 150,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   set.seed(54321)
#   result2 <- SMC_turcotte_cpp(
#     english_words$Yvec,
#     english_words$Tvec,
#     length_UI = 0.5,
#     n_particle = n_particle,
#     U = english_words$U,
#     W = english_words$W,
#     K = english_words$K,
#     num_logs = english_words$num_logs,
#     lambdamat = english_words$lambdamat,
#     keyvec = english_words$keyvec,
#     etavec = english_words$etavec,
#     key0vec = english_words$key0vec,
#     eta0vec = english_words$eta0vec,
#     alphavec = english_words$alphavec,
#     muvec = english_words$muvec,
#     probvec_V = english_words$probvec_V,
#     probvec_Z = english_words$probvec_Z,
#     probvec_Q = english_words$probvec_Q,
#     probvec_F = english_words$probvec_F,
#     P0 = english_words$P0,
#     minimum_n = english_words$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 150,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   # Results should be identical
#   expect_equal(result1$n_UI, result2$n_UI)
#   expect_equal(result1$UI_bounds, result2$UI_bounds)
#   expect_equal(result1$storage_S, result2$storage_S)
  
#   # Check a few particles match
#   for (j in 1:min(3, result1$n_UI)) {
#     expect_equal(result1$storage_S[[j]], result2$storage_S[[j]])
#     for (p in 1:min(3, n_particle)) {
#       expect_equal(result1$storage_B[[j]][[p]], result2$storage_B[[j]][[p]])
#       expect_equal(result1$storage_V[[j]][[p]], result2$storage_V[[j]][[p]])
#     }
#   }
# })

# test_that("SMC_turcotte_cpp handles different Smax values with english_words", {
#   english_words <- load_english_words_data()
  
#   if (is.null(english_words)) {
#     skip("english_words data not available")
#   }
  
#   n_particle <- 50
  
#   for (Smax in c(50, 100, 150)) {
#     result <- SMC_turcotte_cpp(
#       english_words$Yvec,
#       english_words$Tvec,
#       length_UI = 0.5,
#       n_particle = n_particle,
#       U = english_words$U,
#       W = english_words$W,
#       K = english_words$K,
#       num_logs = english_words$num_logs,
#       lambdamat = english_words$lambdamat,
#       keyvec = english_words$keyvec,
#       etavec = english_words$etavec,
#       key0vec = english_words$key0vec,
#       eta0vec = english_words$eta0vec,
#       alphavec = english_words$alphavec,
#       muvec = english_words$muvec,
#       probvec_V = english_words$probvec_V,
#       probvec_Z = english_words$probvec_Z,
#       probvec_Q = english_words$probvec_Q,
#       probvec_F = english_words$probvec_F,
#       P0 = english_words$P0,
#       minimum_n = english_words$minimum_n,
#       Jss1 = 1/3,
#       Js1s = 1/3,
#       Smax = Smax,
#       n_ite = 1000,
#       burn_in = 100,
#       thinning = 5
#     )
    
#     validate_SMC_turcotte_output(result, n_particle)
    
#     # Check that no particle exceeds Smax
#     for (j in 1:result$n_UI) {
#       expect_true(all(result$storage_S[[j]] <= Smax))
#     }
#   }
# })

# test_that("SMC_turcotte_cpp handles subset of english_words data", {
#   english_words <- load_english_words_data()
  
#   if (is.null(english_words)) {
#     skip("english_words data not available")
#   }
  
#   # Use only first 200 observations for faster testing
#   n_subset <- min(200, length(english_words$Yvec))
#   Yvec_subset <- english_words$Yvec[1:n_subset]
#   Tvec_subset <- english_words$Tvec[1:n_subset]
  
#   n_particle <- 30
  
#   result <- SMC_turcotte_cpp(
#     Yvec_subset,
#     Tvec_subset,
#     length_UI = 0.5,
#     n_particle = n_particle,
#     U = english_words$U,
#     W = english_words$W,
#     K = english_words$K,
#     num_logs = english_words$num_logs,
#     lambdamat = english_words$lambdamat,
#     keyvec = english_words$keyvec,
#     etavec = english_words$etavec,
#     key0vec = english_words$key0vec,
#     eta0vec = english_words$eta0vec,
#     alphavec = english_words$alphavec,
#     muvec = english_words$muvec,
#     probvec_V = english_words$probvec_V,
#     probvec_Z = english_words$probvec_Z,
#     probvec_Q = english_words$probvec_Q,
#     probvec_F = english_words$probvec_F,
#     P0 = english_words$P0,
#     minimum_n = english_words$minimum_n,
#     Jss1 = 1/3,
#     Js1s = 1/3,
#     Smax = 100,
#     n_ite = 1000,
#     burn_in = 100,
#     thinning = 5
#   )
  
#   validate_SMC_turcotte_output(result, n_particle)
# })


####### TEST REAL DATA

test_that("SMC_turcotte_cpp returns correct structure with real_data data", {
  
  n_particle <- 2000
  
  result <- SMC_turcotte_cpp(
    real_data$Yvec,
    real_data$Tvec,
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
    n_ite = 20000,
    burn_in = 5000,
    thinning = 5
  )
  
  validate_SMC_turcotte_output(result, n_particle)
})
