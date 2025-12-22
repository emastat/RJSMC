# Unit tests for C++ functions (via R interface)

context("C++ Functions")

test_that("SMC_turcotte_cpp runs without error", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  expect_error({
    result <- SMC_turcotte_cpp(
      ts_data$Yvec, ts_data$Tvec,
      length_UI = 5.0,
      n_particle = 50,
      U = parameters$U,
      W = parameters$W,
      K = parameters$K,
      num_logs = 5,
      lambdamat = parameters$lambdamat,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      P0 = parameters$P0,
      minimum_n = parameters$minimum_n,
      Jss1 = 1/3,
      Js1s = 1/3,
      Smax = 50,
      n_ite = 1000,
      burn_in = 100,
      thinning = 5
    )
  }, NA)
  
  # Check result structure
  result <- SMC_turcotte_cpp(
    ts_data$Yvec, ts_data$Tvec,
    length_UI = 5.0, n_particle = 50,
    U = parameters$U, W = parameters$W, K = parameters$K, num_logs = 5,
    lambdamat = parameters$lambdamat,
    keyvec = parameters$keyvec, etavec = parameters$etavec,
    key0vec = parameters$key0vec, eta0vec = parameters$eta0vec,
    alphavec = parameters$alphavec, muvec = parameters$muvec,
    probvec_V = parameters$probvec_V, probvec_Z = parameters$probvec_Z,
    probvec_Q = parameters$probvec_Q, probvec_F = parameters$probvec_F,
    P0 = parameters$P0, minimum_n = parameters$minimum_n,
    Jss1 = 1/3, Js1s = 1/3, Smax = 50,
    n_ite = 1000, burn_in = 100, thinning = 5
  )
  
  expect_type(result, "list")
  expect_true("n_UI" %in% names(result))
  expect_true("UI_bounds" %in% names(result))
  expect_true("storage_B" %in% names(result))
  expect_true("storage_V" %in% names(result))
  expect_true("storage_weight" %in% names(result))
})

test_that("compute_posterior works correctly", {
  num_discr_intervals <- 10
  num_particles <- 5
  num_states <- 3
  
  # Create state container: each column is a particle, each row is a time point
  state_container <- matrix(
    sample(1:num_states, num_discr_intervals * num_particles, replace = TRUE),
    nrow = num_discr_intervals,
    ncol = num_particles
  )
  
  # Create weights (should be normalized)
  weight_vec <- rep(1.0 / num_particles, num_particles)
  
  result <- compute_posterior(
    num_discr_intervals,
    num_particles,
    state_container,
    num_states,
    weight_vec
  )
  
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), num_discr_intervals)
  expect_equal(ncol(result), num_states)
  
  # Check that rows sum to 1 (probability distribution)
  row_sums <- rowSums(result)
  expect_true(all(abs(row_sums - 1.0) < 1e-10))
})

test_that("extract_index finds correct indices", {
  Tvec <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
  LB <- 1.2
  UB <- 3.2
  
  result <- extract_index(Tvec, LB, UB)
  
  expect_type(result, "double")
  expect_length(result, 3)
  
  # result[1] should be count, result[2] start index, result[3] end index
  expect_true(result[1] >= 0)  # Count
  expect_true(result[2] >= 0)  # Start index
  expect_true(result[3] >= 0)  # End index
  expect_true(result[2] <= result[3])  # Start <= End
})

test_that("cpp_findInterval works like R's findInterval", {
  x <- c(0.5, 1.5, 2.5, 3.5, 4.5)
  breaks <- c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0)
  
  result_cpp <- cpp_findInterval(x, breaks)
  result_r <- findInterval(x, breaks)
  
  expect_equal(result_cpp, result_r)
})

test_that("full_conditional_V samples correctly", {
  T_seg <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  Y_seg <- c(1, 2, 1, 2, 1)
  U <- 3
  num_logs <- 5
  LB <- 0.0
  UB <- 3.0
  lambdamat <- matrix(rep(1/num_logs, U * num_logs), nrow = U, ncol = num_logs)
  probvec_V <- rep(1/U, U)
  V_left <- 1
  P0 <- 0.1
  
  # Test evaluation (sample_V = FALSE)
  result <- full_conditional_V(
    V = 1, U, num_logs, T_seg, Y_seg,
    LB, UB, lambdamat, probvec_V, V_left, P0,
    end_point = 0.0, open_segment = FALSE, sample_V = FALSE
  )
  
  expect_type(result, "list")
  expect_true("V_new" %in% names(result))
  expect_true("eval_densV" %in% names(result))
  expect_true(result$V_new >= 0 && result$V_new <= U)
  expect_true(is.finite(result$eval_densV))
  
  # Test sampling (sample_V = TRUE)
  result2 <- full_conditional_V(
    V = 0, U, num_logs, T_seg, Y_seg,
    LB, UB, lambdamat, probvec_V, V_left, P0,
    end_point = 0.0, open_segment = FALSE, sample_V = TRUE
  )
  
  expect_true(result2$V_new >= 0 && result2$V_new <= U)
})

test_that("full_conditional_Z works correctly", {
  K <- 2
  L <- 2.5
  keyvec <- c(2.0, 3.0)
  etavec <- c(1.0, 1.5)
  probvec_Z <- c(0.5, 0.5)
  
  # Test evaluation
  result <- full_conditional_Z(
    Z = 1, K, L, keyvec, etavec, probvec_Z, sample_Z = FALSE
  )
  
  expect_type(result, "list")
  expect_true("Z_new" %in% names(result))
  expect_true("eval_densZ" %in% names(result))
  expect_true(result$Z_new >= 1 && result$Z_new <= K)
  
  # Test sampling
  result2 <- full_conditional_Z(
    Z = 0, K, L, keyvec, etavec, probvec_Z, sample_Z = TRUE
  )
  
  expect_true(result2$Z_new >= 1 && result2$Z_new <= K)
})

test_that("full_conditional_F works correctly", {
  F_val <- 1
  L <- 2.0
  key0vec <- c(2.0, 3.0)
  eta0vec <- c(1.0, 1.5)
  probvec_F <- c(0.5, 0.5)
  
  result <- full_conditional_F(
    F_val, L, key0vec, eta0vec, probvec_F, sample_F = FALSE
  )
  
  expect_type(result, "list")
  expect_true("F_new" %in% names(result))
  expect_true("eval_densF" %in% names(result))
  expect_true(result$F_new >= 1 && result$F_new <= 2)
})

test_that("incremental_weight computes correctly", {
  T_seg <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
  Y_seg <- c(1, 2, 1, 2, 1, 2)
  V <- 1
  Z <- 1
  Q <- 1
  F_val <- 0
  U <- 3
  K <- 2
  W <- 2
  B_last <- 0.0
  B1 <- 2.0
  num_logs <- 5
  lambdamat <- matrix(rep(1/num_logs, U * num_logs), nrow = U, ncol = num_logs)
  keyvec <- c(2.0, 3.0)
  etavec <- c(1.0, 1.5)
  key0vec <- c(2.0, 3.0)
  eta0vec <- c(1.0, 1.5)
  alphavec <- c(2.0, 3.0)
  muvec <- c(1.0, 1.5)
  probvec_V <- rep(1/U, U)
  probvec_Z <- rep(1/K, K)
  probvec_Q <- rep(1/W, W)
  probvec_F <- c(0.5, 0.5)
  V_last_complete <- 1
  P0 <- 0.1
  minimum_n <- 0
  t_star <- 0.5
  
  result <- incremental_weight(
    T_seg, Y_seg, V, Z, Q, F_val, U, K, W,
    B_last, B1, num_logs, lambdamat,
    keyvec, etavec, key0vec, eta0vec,
    alphavec, muvec,
    probvec_V, probvec_Z, probvec_Q, probvec_F,
    V_last_complete, P0, minimum_n, t_star
  )
  
  expect_type(result, "double")
  expect_true(is.finite(result))
  # Incremental weight should be a log value (can be negative)
})

