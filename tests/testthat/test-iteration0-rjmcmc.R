# Tests for iteration0_RJMCMC function

context("iteration0_RJMCMC")

test_that("iteration0_RJMCMC runs without error with valid inputs", {
  # Create test data
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 123)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  # Extract segment data within an interval
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  # Filter observations within the interval
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  # Only proceed if we have observations
  if (length(T_seg) > 0) {
    set.seed(456)
    result <- iteration0_RJMCMC(
      T_seg = T_seg,
      Y_seg = Y_seg,
      minimum_n = parameters$minimum_n,
      start_point = start_point,
      end_point = end_point,
      t_star = t_star,
      K = parameters$K,
      W = parameters$W,
      U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    expect_type(result, "list")
  }
})

test_that("iteration0_RJMCMC returns correct structure", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 789)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    set.seed(111)
    result <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Check all required elements exist
    expect_true("S" %in% names(result))
    expect_true("Bvec" %in% names(result))
    expect_true("Vvec" %in% names(result))
    expect_true("Zvec" %in% names(result))
    expect_true("Qvec" %in% names(result))
    expect_true("Fvec" %in% names(result))
    expect_true("Nvec" %in% names(result))
    
    # Check types
    expect_type(result$S, "integer")
    expect_type(result$Bvec, "double")
    expect_type(result$Vvec, "integer")
    expect_type(result$Zvec, "integer")
    expect_type(result$Qvec, "integer")
    expect_type(result$Fvec, "integer")
    expect_type(result$Nvec, "integer")
  }
})

test_that("iteration0_RJMCMC returns correct dimensions", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 222)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  Smax <- 50
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    set.seed(333)
    result <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = Smax
    )
    
    S <- result$S
    
    # Check dimensions match Smax
    expect_equal(length(result$Bvec), Smax + 1)
    expect_equal(length(result$Vvec), Smax)
    expect_equal(length(result$Zvec), Smax)
    expect_equal(length(result$Qvec), Smax)
    expect_equal(length(result$Fvec), Smax)
    
    # Nvec should have length S (not Smax)
    expect_equal(length(result$Nvec), S)
    
    # S should be positive and <= Smax
    expect_true(S > 0)
    expect_true(S <= Smax)
  }
})

test_that("iteration0_RJMCMC breakpoints are valid", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 444)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    set.seed(555)
    result <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    S <- result$S
    Bvec <- result$Bvec[1:(S+1)]  # Only use the first S+1 breakpoints
    
    # Breakpoints should be sorted
    expect_true(all(diff(Bvec) >= 0))
    
    # First breakpoint should be >= t_star
    expect_true(Bvec[1] >= t_star)
    
    # Last breakpoint should be > end_point (as per documentation)
    expect_true(tail(Bvec, 1) > end_point)
    
    # All breakpoints should be finite
    expect_true(all(is.finite(Bvec)))
  }
})

test_that("iteration0_RJMCMC states are within valid ranges", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 666)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    set.seed(777)
    result <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    S <- result$S
    
    # Extract only the first S elements (rest are padding)
    Vvec <- result$Vvec[1:S]
    Zvec <- result$Zvec[1:S]
    Qvec <- result$Qvec[1:S]
    Fvec <- result$Fvec[1:S]
    
    # V states: 0 for empty segments, 1 to U (1-indexed) for non-empty segments
    # Valid range: [0, U]
    expect_true(all(Vvec >= 0))
    expect_true(all(Vvec <= parameters$U))
    
    # Z states: 0 for empty segments, 1 to K (1-indexed) for non-empty segments
    # Valid range: [0, K]
    expect_true(all(Zvec >= 0))
    expect_true(all(Zvec <= parameters$K))
    
    # Q states: 0 for empty segments, 1 to W (1-indexed) for non-empty segments
    # Valid range: [0, W]
    expect_true(all(Qvec >= 0))
    expect_true(all(Qvec <= parameters$W))
    
    # F states: 0 for non-empty segments, 1 to F_levels (1-indexed) for empty segments
    # Valid range: [0, 2] (typically 2 F states)
    expect_true(all(Fvec >= 0))
    expect_true(all(Fvec <= 2))
    
    # Nvec should be non-negative
    expect_true(all(result$Nvec >= 0))
  }
})

test_that("iteration0_RJMCMC handles empty_mix parameter", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 888)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    # Test with empty_mix = TRUE
    set.seed(999)
    result1 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Test with empty_mix = FALSE
    set.seed(1000)
    result2 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = FALSE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Both should produce valid results
    expect_type(result1, "list")
    expect_type(result2, "list")
    expect_true(result1$S > 0)
    expect_true(result2$S > 0)
  }
})

test_that("iteration0_RJMCMC respects Smax constraint", {
  ts_data <- create_test_ts_data(n_obs = 100, num_logs = 5, seed = 1111)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 20.0  # Larger interval to potentially generate more segments
  t_star <- 0.0
  Smax <- 10  # Small Smax to test constraint
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    set.seed(2222)
    result <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = Smax
    )
    
    # S should not exceed Smax
    expect_true(result$S <= Smax)
    
    # Vectors should still have Smax dimensions
    expect_equal(length(result$Vvec), Smax)
    expect_equal(length(result$Zvec), Smax)
    expect_equal(length(result$Qvec), Smax)
    expect_equal(length(result$Fvec), Smax)
    expect_equal(length(result$Bvec), Smax + 1)
  }
})

test_that("iteration0_RJMCMC handles minimum_n constraint", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 3333)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    # Test with minimum_n = 0 (no constraint)
    set.seed(4444)
    result1 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = 0,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Test with minimum_n = 3 (more restrictive)
    set.seed(5555)
    result2 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = 3,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Both should produce valid results
    expect_type(result1, "list")
    expect_type(result2, "list")
    
    # With higher minimum_n, we might get fewer segments
    # (but this depends on the data, so we just check validity)
    expect_true(result1$S > 0)
    expect_true(result2$S > 0)
  }
})

test_that("iteration0_RJMCMC produces consistent results with same seed", {
  ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 6666)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 10.0
  t_star <- 0.0
  
  in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
  T_seg <- ts_data$Tvec[in_interval]
  Y_seg <- ts_data$Yvec[in_interval]
  
  if (length(T_seg) > 0) {
    # Run twice with same seed
    set.seed(7777)
    result1 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    set.seed(7777)
    result2 <- iteration0_RJMCMC(
      T_seg, Y_seg,
      minimum_n = parameters$minimum_n,
      start_point, end_point, t_star,
      K = parameters$K, W = parameters$W, U = parameters$U,
      empty_mix = TRUE,
      probvec_V = parameters$probvec_V,
      probvec_Z = parameters$probvec_Z,
      probvec_Q = parameters$probvec_Q,
      probvec_F = parameters$probvec_F,
      alphavec = parameters$alphavec,
      muvec = parameters$muvec,
      keyvec = parameters$keyvec,
      etavec = parameters$etavec,
      key0vec = parameters$key0vec,
      eta0vec = parameters$eta0vec,
      lambdamat = parameters$lambdamat,
      P0 = parameters$P0,
      num_logs = 5,
      max_range = 1.0,
      Smax = 50
    )
    
    # Results should be identical with same seed
    expect_equal(result1$S, result2$S)
    expect_equal(result1$Bvec, result2$Bvec)
    expect_equal(result1$Vvec, result2$Vvec)
    expect_equal(result1$Zvec, result2$Zvec)
    expect_equal(result1$Qvec, result2$Qvec)
    expect_equal(result1$Fvec, result2$Fvec)
    expect_equal(result1$Nvec, result2$Nvec)
  }
})

test_that("iteration0_RJMCMC handles edge case with few observations", {
  # Create data with very few observations
  T_seg <- c(1.0, 2.0, 3.0)
  Y_seg <- c(1, 2, 1)
  parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
  
  start_point <- 0.0
  end_point <- 5.0
  t_star <- 0.0
  
  set.seed(8888)
  result <- iteration0_RJMCMC(
    T_seg, Y_seg,
    minimum_n = 0,
    start_point, end_point, t_star,
    K = parameters$K, W = parameters$W, U = parameters$U,
    empty_mix = TRUE,
    probvec_V = parameters$probvec_V,
    probvec_Z = parameters$probvec_Z,
    probvec_Q = parameters$probvec_Q,
    probvec_F = parameters$probvec_F,
    alphavec = parameters$alphavec,
    muvec = parameters$muvec,
    keyvec = parameters$keyvec,
    etavec = parameters$etavec,
    key0vec = parameters$key0vec,
    eta0vec = parameters$eta0vec,
    lambdamat = parameters$lambdamat,
    P0 = parameters$P0,
    num_logs = 5,
    max_range = 1.0,
    Smax = 50
  )
  
  # Should still produce valid output
  expect_type(result, "list")
  expect_true(result$S > 0)
  expect_equal(length(result$Nvec), result$S)
})

