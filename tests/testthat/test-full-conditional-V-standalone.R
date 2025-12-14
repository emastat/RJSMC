# Standalone test for full_conditional_V to identify out-of-bounds access
# This tests the function with T_seg and Y_seg of size 16 (matching the error case)

context("full_conditional_V - Standalone Debug Test")

test_that("full_conditional_V with T_seg and Y_seg size 16", {
  # Test Case: T_seg and Y_seg have 16 elements
  # This matches the scenario where we get "index 16 >= vector size 16"
  set.seed(123)
  T_seg <- cumsum(rexp(16, rate = 1))  # 16 time points
  Y_seg <- sample(1:5, size = 16, replace = TRUE)  # 16 messages
  
  # Create parameters
  U <- 3  # Number of V states
  num_logs <- 5  # Number of message types
  lambdamat <- matrix(runif(U * num_logs), nrow = U, ncol = num_logs)
  lambdamat <- lambdamat / rowSums(lambdamat)  # Normalize rows
  probvec_V <- c(0.3, 0.4, 0.3)  # Probability for each V state
  P0 <- 0.1
  V_left <- 1
  
  cat("\n=== Test 1: T_seg and Y_seg size 16 ===\n")
  cat("T_seg size:", length(T_seg), "\n")
  cat("Y_seg size:", length(Y_seg), "\n")
  cat("T_seg range:", min(T_seg), "to", max(T_seg), "\n")
  
  # Test with a segment that covers all observations
  LB <- min(T_seg) - 0.1
  UB <- max(T_seg) + 0.1
  
  cat("LB:", LB, ", UB:", UB, "\n")
  
  tryCatch({
    result <- full_conditional_V(
      V = 0,
      U = U,
      num_logs = num_logs,
      T_seg = T_seg,
      Y_seg = Y_seg,
      LB = LB,
      UB = UB,
      lambdamat = lambdamat,
      probvec_V = probvec_V,
      V_left = V_left,
      P0 = P0,
      end_point = 0.0,
      open_segment = FALSE,
      sample_V = FALSE
    )
    
    cat("Result V:", result$V, "\n")
    cat("Result eval_densV:", result$eval_densV, "\n")
    
    expect_true("V" %in% names(result))
    expect_true("eval_densV" %in% names(result))
    
    cat("✓ Test 1 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 1 FAILED:", e$message, "\n")
    stop(e)
  })
})

test_that("full_conditional_V with segment at end of T_seg", {
  # Test Case: Segment that includes the last observation
  # This could trigger index 16 if extract_index returns end_idx = 16
  set.seed(456)
  T_seg2 <- cumsum(rexp(16, rate = 1))
  Y_seg2 <- sample(1:5, size = 16, replace = TRUE)
  
  U <- 3
  num_logs <- 5
  lambdamat2 <- matrix(runif(U * num_logs), nrow = U, ncol = num_logs)
  lambdamat2 <- lambdamat2 / rowSums(lambdamat2)
  probvec_V2 <- c(0.3, 0.4, 0.3)
  P0 <- 0.1
  V_left <- 1
  
  cat("\n=== Test 2: Segment at end of T_seg ===\n")
  cat("T_seg2 size:", length(T_seg2), "\n")
  cat("Y_seg2 size:", length(Y_seg2), "\n")
  cat("T_seg2 range:", min(T_seg2), "to", max(T_seg2), "\n")
  
  # Segment covering last few observations
  LB <- T_seg2[14] - 0.01  # Just before second-to-last observation
  UB <- T_seg2[15] + 0.01  # Just after last observation (index 15)
  
  cat("LB:", LB, ", UB:", UB, "\n")
  cat("Last T_seg2 value:", T_seg2[15], "\n")
  
  tryCatch({
    result2 <- full_conditional_V(
      V = 0,
      U = U,
      num_logs = num_logs,
      T_seg = T_seg2,
      Y_seg = Y_seg2,
      LB = LB,
      UB = UB,
      lambdamat = lambdamat2,
      probvec_V = probvec_V2,
      V_left = V_left,
      P0 = P0,
      end_point = 0.0,
      open_segment = FALSE,
      sample_V = FALSE
    )
    
    cat("Result2 V:", result2$V, "\n")
    expect_true("V" %in% names(result2))
    
    cat("✓ Test 2 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 2 FAILED:", e$message, "\n")
    stop(e)
  })
})

test_that("full_conditional_V with open_segment=true", {
  # Test Case: Open segment (partially observed)
  set.seed(789)
  T_seg3 <- cumsum(rexp(16, rate = 1))
  Y_seg3 <- sample(1:5, size = 16, replace = TRUE)
  
  U <- 3
  num_logs <- 5
  lambdamat3 <- matrix(runif(U * num_logs), nrow = U, ncol = num_logs)
  lambdamat3 <- lambdamat3 / rowSums(lambdamat3)
  probvec_V3 <- c(0.3, 0.4, 0.3)
  P0 <- 0.1
  V_left <- 1
  
  cat("\n=== Test 3: Open segment (partially observed) ===\n")
  cat("T_seg3 size:", length(T_seg3), "\n")
  cat("Y_seg3 size:", length(Y_seg3), "\n")
  
  LB <- T_seg3[10] - 0.01
  end_point <- T_seg3[15] + 0.01  # end_point beyond last observation
  
  cat("LB:", LB, ", end_point:", end_point, "\n")
  
  tryCatch({
    result3 <- full_conditional_V(
      V = 0,
      U = U,
      num_logs = num_logs,
      T_seg = T_seg3,
      Y_seg = Y_seg3,
      LB = LB,
      UB = 0.0,  # Not used when open_segment=true
      lambdamat = lambdamat3,
      probvec_V = probvec_V3,
      V_left = V_left,
      P0 = P0,
      end_point = end_point,
      open_segment = TRUE,
      sample_V = FALSE
    )
    
    cat("Result3 V:", result3$V, "\n")
    expect_true("V" %in% names(result3))
    
    cat("✓ Test 3 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 3 FAILED:", e$message, "\n")
    stop(e)
  })
})

test_that("full_conditional_V exact scenario from iteration0_RJMCMC", {
  # Test Case: Exact scenario matching the error case
  # This simulates what happens when initial_state calls full_conditional_V
  set.seed(123)
  T_seg4 <- cumsum(rexp(16, rate = 1))
  Y_seg4 <- sample(1:5, size = 16, replace = TRUE)
  
  # Simulate Bvec with 16 breakpoints (S=15 segments)
  # The last segment would use Bvec[14] and Bvec[15] as LB and UB
  Bvec4 <- c(0.0, T_seg4[1:14], max(T_seg4) + 1.0)
  
  U <- 3
  num_logs <- 5
  lambdamat4 <- matrix(runif(U * num_logs), nrow = U, ncol = num_logs)
  lambdamat4 <- lambdamat4 / rowSums(lambdamat4)
  probvec_V4 <- c(0.3, 0.4, 0.3)
  P0 <- 0.1
  V_left <- 1
  
  cat("\n=== Test 4: Exact scenario from iteration0_RJMCMC ===\n")
  cat("T_seg4 size:", length(T_seg4), "\n")
  cat("Y_seg4 size:", length(Y_seg4), "\n")
  cat("Bvec4 size:", length(Bvec4), "\n")
  cat("S (segments):", length(Bvec4) - 1, "\n")
  
  # Test the last segment (i=14, uses Bvec[14] and Bvec[15])
  i <- 14
  LB <- Bvec4[i]
  UB <- Bvec4[i+1]
  
  cat("Testing segment i=", i, "\n")
  cat("LB (Bvec[", i, "]):", LB, "\n")
  cat("UB (Bvec[", i+1, "]):", UB, "\n")
  
  tryCatch({
    result4 <- full_conditional_V(
      V = 0,
      U = U,
      num_logs = num_logs,
      T_seg = T_seg4,
      Y_seg = Y_seg4,
      LB = LB,
      UB = UB,
      lambdamat = lambdamat4,
      probvec_V = probvec_V4,
      V_left = V_left,
      P0 = P0,
      end_point = 0.0,
      open_segment = FALSE,
      sample_V = FALSE
    )
    
    cat("Result4 V:", result4$V, "\n")
    expect_true("V" %in% names(result4))
    
    cat("✓ Test 4 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 4 FAILED:", e$message, "\n")
    cat("  This is the exact scenario causing the out-of-bounds warning!\n")
    stop(e)
  })
})

