# Tests for breakpoints_sampling function

context("breakpoints_sampling")

library(mclust)
library(fitdistrplus)

test_that("breakpoints_sampling works correctly with valid inputs", {
  # Create a list of breakpoint vectors (simulating RJMCMC output)
  start_point <- 0.0
  end_point <- 10.0
  sample_size <- 5
  
  # Create breakpoint list with multiple particles
  breakpoint_list <- list(
    c(2.0, 5.0, 8.0),
    c(1.5, 4.5, 7.5),
    c(2.5, 6.0, 9.0),
    c(1.0, 3.0, 5.0, 7.0),
    c(2.0, 4.0, 6.0, 8.0)
  )
  
  set.seed(123)  # For reproducibility
  result <- breakpoints_sampling(
    start_point = start_point,
    end_point = end_point,
    breakpoint_list = breakpoint_list,
    sample_size = sample_size
  )
  
  # Check output structure
  expect_type(result, "list")
  expect_true("log_density" %in% names(result))
  expect_true("S.samp" %in% names(result))
  
  # Check log_density
  expect_type(result$log_density, "double")
  expect_length(result$log_density, sample_size)
  expect_false(any(is.na(result$log_density)))
  expect_false(any(is.nan(result$log_density)))
  
  # Check S.samp
  expect_type(result$S.samp, "list")
  expect_length(result$S.samp, sample_size)
  
  # Check each sampled breakpoint vector
  for (i in 1:sample_size) {
    expect_type(result$S.samp[[i]], "double")
    expect_false(any(is.na(result$S.samp[[i]])))
    expect_false(any(is.nan(result$S.samp[[i]])))
    
    # Check breakpoints are within bounds
    expect_true(all(result$S.samp[[i]] >= start_point))
    expect_true(all(result$S.samp[[i]] <= end_point))
    
    # Check breakpoints are sorted
    expect_true(all(diff(result$S.samp[[i]]) >= 0))
    
    # Check that start_point and end_point are included
    expect_true(result$S.samp[[i]][1] == start_point)
    expect_true(tail(result$S.samp[[i]], 1) == end_point)
  }
})

test_that("breakpoints_sampling handles single breakpoint case", {
  start_point <- 0.0
  end_point <- 10.0
  sample_size <- 3
  
  # Multiple particles with the same single breakpoint
  # This is more realistic than a single particle and avoids Mclust hanging
  # with too few data points
  breakpoint_list <- list(
    c(5.0),
    c(5.0),
    c(5.0),
    c(5.0),
    c(5.0)
  )
  
  set.seed(456)
  result <- breakpoints_sampling(
    start_point = start_point,
    end_point = end_point,
    breakpoint_list = breakpoint_list,
    sample_size = sample_size
  )
  
  expect_type(result, "list")
  expect_length(result$log_density, sample_size)
  expect_length(result$S.samp, sample_size)
  
  # Each sample should have at least start and end points
  for (i in 1:sample_size) {
    expect_true(length(result$S.samp[[i]]) >= 2)
    expect_true(result$S.samp[[i]][1] == start_point)
    expect_true(tail(result$S.samp[[i]], 1) == end_point)
  }
})

test_that("breakpoints_sampling handles multiple breakpoints per particle", {
  start_point <- 0.0
  end_point <- 20.0
  sample_size <- 4
  
  # Breakpoint list with varying numbers of breakpoints
  breakpoint_list <- list(
    c(5.0, 10.0),
    c(3.0, 7.0, 12.0, 15.0),
    c(2.0, 4.0, 6.0, 8.0, 10.0, 12.0),
    c(1.0, 5.0, 9.0)
  )
  
  set.seed(789)
  result <- breakpoints_sampling(
    start_point = start_point,
    end_point = end_point,
    breakpoint_list = breakpoint_list,
    sample_size = sample_size
  )
  
  expect_length(result$S.samp, sample_size)
  expect_length(result$log_density, sample_size)
  
  # Verify all breakpoints are valid
  for (i in 1:sample_size) {
    expect_true(all(result$S.samp[[i]] >= start_point))
    expect_true(all(result$S.samp[[i]] <= end_point))
    expect_true(all(diff(result$S.samp[[i]]) >= 0))
  }
})

test_that("breakpoints_sampling filters breakpoints outside interval", {
  start_point <- 5.0
  end_point <- 15.0
  sample_size <- 3
  
  # Breakpoint list with some breakpoints outside the interval
  breakpoint_list <- list(
    c(2.0, 7.0, 10.0, 18.0),  # 2.0 and 18.0 are outside
    c(6.0, 8.0, 12.0),
    c(5.5, 9.0, 14.5)
  )
  
  set.seed(321)
  result <- breakpoints_sampling(
    start_point = start_point,
    end_point = end_point,
    breakpoint_list = breakpoint_list,
    sample_size = sample_size
  )
  
  # All sampled breakpoints should be within [start_point, end_point]
  for (i in 1:sample_size) {
    expect_true(all(result$S.samp[[i]] >= start_point))
    expect_true(all(result$S.samp[[i]] <= end_point))
  }
})

test_that("breakpoints_sampling handles edge case with minimal breakpoints", {
  start_point <- 0.0
  end_point <- 10.0
  sample_size <- 2
  
  # Very few breakpoints (edge case)
  breakpoint_list <- list(
    c(5.0),
    c(5.0)  # Only one unique breakpoint
  )
  
  set.seed(654)
  result <- breakpoints_sampling(
    start_point = start_point,
    end_point = end_point,
    breakpoint_list = breakpoint_list,
    sample_size = sample_size
  )
  
  expect_length(result$S.samp, sample_size)
  expect_length(result$log_density, sample_size)
  
  # Should still produce valid output
  for (i in 1:sample_size) {
    expect_true(length(result$S.samp[[i]]) >= 2)
    expect_true(result$S.samp[[i]][1] == start_point)
    expect_true(tail(result$S.samp[[i]], 1) == end_point)
  }
})

test_that("breakpoints_sampling produces different samples", {
  start_point <- 0.0
  end_point <- 10.0
  sample_size <- 10
  
  breakpoint_list <- list(
    c(2.0, 5.0, 8.0),
    c(1.5, 4.5, 7.5),
    c(2.5, 6.0, 9.0)
  )
  
  set.seed(111)
  result1 <- breakpoints_sampling(
    start_point, end_point, breakpoint_list, sample_size
  )
  
  set.seed(222)
  result2 <- breakpoints_sampling(
    start_point, end_point, breakpoint_list, sample_size
  )
  
  # Results should be different with different seeds
  # (at least some breakpoints should differ)
  all_same <- TRUE
  for (i in 1:sample_size) {
    if (!identical(result1$S.samp[[i]], result2$S.samp[[i]])) {
      all_same <- FALSE
      break
    }
  }
  expect_false(all_same)  # Should produce different samples
})

test_that("breakpoints_sampling log_density values are reasonable", {
  start_point <- 0.0
  end_point <- 10.0
  sample_size <- 5
  
  breakpoint_list <- list(
    c(2.0, 5.0, 8.0),
    c(1.5, 4.5, 7.5),
    c(2.5, 6.0, 9.0)
  )
  
  set.seed(999)
  result <- breakpoints_sampling(
    start_point, end_point, breakpoint_list, sample_size
  )
  
  # Log densities should be finite and not too extreme
  expect_true(all(is.finite(result$log_density)))
  expect_true(all(result$log_density <= 0))  # Log probabilities are <= 0
  expect_true(all(result$log_density >= -100))  # Function caps at -100
})

