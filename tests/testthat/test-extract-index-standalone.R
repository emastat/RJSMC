# Standalone test for extract_index to identify out-of-bounds access
# This tests the function with x of size 16 (matching the error case)

context("extract_index - Standalone Debug Test")

test_that("extract_index with x size 16 - Test Case 1", {
  # Test Case: x has 16 elements
  # This matches the scenario where we get "index 16 >= vector size 16"
  set.seed(123)
  x <- cumsum(rexp(16, rate = 1))  # 16 elements
  
  cat("\n=== Test 1: x size 16 ===\n")
  cat("x size:", length(x), "\n")
  cat("x range:", min(x), "to", max(x), "\n")
  
  # Test with bounds that include all elements
  lower <- min(x) - 0.1
  upper <- max(x) + 0.1
  
  cat("lower:", lower, ", upper:", upper, "\n")
  
  result <- extract_index(x, lower, upper)
  
  # NOTE: R uses 1-based indexing, C++ uses 0-based
  # C++ out[0] = length, out[1] = start index, out[2] = end index
  # R result[1] = length, result[2] = start index, result[3] = end index
  cat("Result length (out[0] / result[1]):", result[1], "\n")
  cat("Result start index (out[1] / result[2]):", result[2], "\n")
  cat("Result end index (out[2] / result[3]):", result[3], "\n")
  
  # Check bounds - result[2] is the start index (C++ out[1])
  # result[3] is the end index (C++ out[2])
  expect_true(result[3] >= -1)  # -1 means no elements found
  expect_true(result[3] < length(x))  # Must be < x.size() for valid index
  expect_true(result[2] >= 0)  # Start index (C++ out[1])
  expect_true(result[2] < length(x))  # Start index must be < x.size()
  
  if(result[2] >= length(x)){
    cat("✗ PROBLEM: start index (result[2]) =", result[2], ">= x.size() =", length(x), "\n")
    stop("extract_index returns out-of-bounds start index!")
  }
  if(result[3] >= length(x)){
    cat("✗ PROBLEM: end index (result[3]) =", result[3], ">= x.size() =", length(x), "\n")
    stop("extract_index returns out-of-bounds end index!")
  }
  
  cat("✓ All indices are valid\n")
  
  cat("✓ Test 1 PASSED\n")
})

test_that("extract_index with bounds at end of vector", {
  # Test Case: Upper bound just after last element
  # This could trigger index 16 if not handled correctly
  set.seed(456)
  x2 <- cumsum(rexp(16, rate = 1))
  
  cat("\n=== Test 2: Bounds at end of vector ===\n")
  cat("x2 size:", length(x2), "\n")
  cat("Last x2 value:", x2[15], "\n")
  
  lower2 <- x2[14] - 0.01  # Just before second-to-last
  upper2 <- x2[15] + 0.01  # Just after last element (index 15)
  
  cat("lower2:", lower2, ", upper2:", upper2, "\n")
  
  result2 <- extract_index(x2, lower2, upper2)
  
  # NOTE: result2[2] is start index, result2[3] is end index
  cat("Result2 start index (result2[2]):", result2[2], "\n")
  cat("Result2 end index (result2[3]):", result2[3], "\n")
  cat("x2.size():", length(x2), "\n")
  
  # Critical check: end index (result2[3]) must be < x2.size()
  if(result2[3] >= length(x2)){
    cat("✗ PROBLEM: end index (result2[3]) =", result2[3], ">= x2.size() =", length(x2), "\n")
    stop("extract_index returns out-of-bounds index at end of vector!")
  } else {
    cat("✓ End index is valid\n")
  }
  
  expect_true(result2[3] < length(x2))
  
  cat("✓ Test 2 PASSED\n")
})

test_that("extract_index with upper bound beyond last element", {
  # Test Case: Upper bound way beyond last element
  set.seed(789)
  x3 <- cumsum(rexp(16, rate = 1))
  
  cat("\n=== Test 3: Upper bound beyond last element ===\n")
  cat("x3 size:", length(x3), "\n")
  cat("Last x3 value:", x3[15], "\n")
  
  # FIX: R uses 1-based indexing, so x3[1] is the first element, not x3[0]
  lower3 <- x3[1] - 0.1  # or could use min(x3) - 0.1
  upper3 <- x3[15] + 100.0  # Way beyond last element
  
  cat("lower3:", lower3, ", upper3:", upper3, "\n")
  
  result3 <- extract_index(x3, lower3, upper3)
  
  # NOTE: result3[2] is start index, result3[3] is end index
  cat("Result3 start index (result3[2]):", result3[2], "\n")
  cat("Result3 end index (result3[3]):", result3[3], "\n")
  cat("x3.size():", length(x3), "\n")
  
  # Should return index 15 (last valid index), not 16
  if(result3[3] >= length(x3)){
    cat("✗ PROBLEM: end index (result3[3]) =", result3[3], ">= x3.size() =", length(x3), "\n")
    stop("extract_index returns out-of-bounds index with large upper bound!")
  } else {
    cat("✓ End index is valid (should be", length(x3) - 1, ")\n")
  }
  
  expect_true(result3[3] < length(x3))
  expect_equal(result3[3], length(x3) - 1)  # Should be last valid index
  
  cat("✓ Test 3 PASSED\n")
})

test_that("extract_index exact scenario from iteration0_RJMCMC", {
  # Test Case: Exact scenario matching the error case
  # T_seg has 16 elements, segment bounds that could cause issue
  set.seed(123)
  T_seg <- cumsum(rexp(16, rate = 1))
  
  # Simulate a segment that includes the last observation
  # Bvec would have 16 breakpoints (S=15 segments)
  # Last segment uses Bvec[14] and Bvec[15] as LB and UB
  Bvec <- c(0.0, T_seg[1:14], max(T_seg) + 1.0)
  
  cat("\n=== Test 4: Exact scenario from iteration0_RJMCMC ===\n")
  cat("T_seg size:", length(T_seg), "\n")
  cat("Bvec size:", length(Bvec), "\n")
  cat("S (segments):", length(Bvec) - 1, "\n")
  
  # Test the last segment (i=14, uses Bvec[14] and Bvec[15])
  i <- 14
  LB <- Bvec[i]
  UB <- Bvec[i+1]
  
  cat("Testing segment i=", i, "\n")
  cat("LB (Bvec[", i, "]):", LB, "\n")
  cat("UB (Bvec[", i+1, "]):", UB, "\n")
  cat("T_seg range:", min(T_seg), "to", max(T_seg), "\n")
  
  result4 <- extract_index(T_seg, LB, UB)
  
  # NOTE: result4[2] is start index, result4[3] is end index
  cat("Result4 start index (result4[2]):", result4[2], "\n")
  cat("Result4 end index (result4[3]):", result4[3], "\n")
  cat("T_seg.size():", length(T_seg), "\n")
  
  # CRITICAL CHECK: end index (result4[3]) must be < T_seg.size()
  if(result4[3] >= length(T_seg)){
    cat("✗ PROBLEM FOUND: end index (result4[3]) =", result4[3], ">= T_seg.size() =", length(T_seg), "\n")
    cat("  This is the exact scenario causing the out-of-bounds warning!\n")
    stop("extract_index returns out-of-bounds index - this is the bug!")
  } else {
    cat("✓ End index is valid\n")
  }
  
  expect_true(result4[3] < length(T_seg))
  
  cat("✓ Test 4 PASSED\n")
})

