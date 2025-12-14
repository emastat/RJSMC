# Standalone test for cpp_findInterval to identify out-of-bounds access
# This tests the function with T_seg of size 16 (matching the error case)

context("cpp_findInterval - Standalone Debug Test")

test_that("cpp_findInterval with T_seg size 16 - Test Case 1", {
  # Test Case 1: T_seg has 16 elements, Bvec has various sizes
  # cpp_findInterval returns 0-based indices in range [0, breaks.size()]
  T_seg <- seq(0.5, 8.0, length.out = 16)  # 16 elements
  Bvec <- c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0)  # 6 elements (S+1 breakpoints, S=5 segments)
  S <- length(Bvec) - 1  # S = 5
  
  cat("\n=== Test 1: T_seg size 16, Bvec size 6 (S=5) ===\n")
  cat("T_seg.size():", length(T_seg), "\n")
  cat("Bvec.size():", length(Bvec), "\n")
  cat("S (number of segments):", S, "\n")
  
  result <- cpp_findInterval(T_seg, Bvec)
  
  cat("Result length:", length(result), "\n")
  cat("Min result value:", min(result), "\n")
  cat("Max result value:", max(result), "\n")
  cat("Unique result values:", sort(unique(result)), "\n")
  
  # cpp_findInterval can return values from 0 to breaks.size() (0 to 6 in this case)
  # This is the valid range for cpp_findInterval itself
  expect_true(all(result >= 0))
  expect_true(all(result <= length(Bvec)))
  
  # CRITICAL CHECK: When used with make_table(interval, S, false), 
  # make_table expects values in [0, S] (0 to 5), but cpp_findInterval can return S+1 (6)
  # This would cause out-of-bounds when make_table tries to access freque[6]
  if(any(result > S)){
    cat("✗ PROBLEM FOUND: cpp_findInterval returned values > S!\n")
    cat("  Invalid values (for make_table):", result[result > S], "\n")
    cat("  S =", S, ", but got values up to", max(result), "\n")
    cat("  This will cause out-of-bounds in make_table(interval, S, false)!\n")
  } else {
    cat("✓ All values are <= S, safe for make_table\n")
  }
  
  cat("✓ Test 1 completed (checking for values > S)\n")
})

test_that("cpp_findInterval - Test Case 2: T_seg values > max(Bvec)", {
  # Edge case - T_seg values outside Bvec range
  # This is the critical case: when T_seg contains values > max(Bvec),
  # cpp_findInterval will return breaks.size() (the maximum valid index)
  T_seg2 <- c(seq(0.1, 4.9, length.out = 15), 10.0)  # Last element > max(Bvec)
  Bvec2 <- c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0)  # 6 elements (S+1 breakpoints, S=5 segments)
  S2 <- length(Bvec2) - 1  # S = 5
  
  cat("\n=== Test 2: T_seg with values > max(Bvec) ===\n")
  cat("T_seg2 size:", length(T_seg2), "\n")
  cat("Bvec2 size:", length(Bvec2), "\n")
  cat("S (number of segments):", S2, "\n")
  cat("T_seg2 range:", min(T_seg2), "to", max(T_seg2), "\n")
  cat("Bvec2 range:", min(Bvec2), "to", max(Bvec2), "\n")
  
  result2 <- cpp_findInterval(T_seg2, Bvec2)
  
  cat("Min result2 value:", min(result2), "\n")
  cat("Max result2 value:", max(result2), "\n")
  cat("Unique result2 values:", sort(unique(result2)), "\n")
  
  # cpp_findInterval can return values from 0 to breaks.size() (0 to 6 in this case)
  # This is the valid range for cpp_findInterval itself
  expect_true(all(result2 >= 0))
  expect_true(all(result2 <= length(Bvec2)))
  
  # CRITICAL CHECK: When T_seg contains values > max(Bvec), cpp_findInterval returns breaks.size()
  # This is expected behavior for cpp_findInterval, but problematic for make_table
  if(any(result2 > S2)){
    cat("✗ PROBLEM FOUND: cpp_findInterval returned values > S!\n")
    cat("  Invalid values (for make_table):", result2[result2 > S2], "\n")
    cat("  S =", S2, ", but got values up to", max(result2), "\n")
    cat("  This will cause out-of-bounds in make_table(interval, S, false)!\n")
    cat("  This is the expected behavior when T_seg contains values > max(Bvec).\n")
  } else {
    cat("✓ All values are <= S, safe for make_table\n")
  }
  
  cat("✓ Test 2 completed (checking for values > S)\n")
})

test_that("cpp_findInterval - Test Case 3: Exact iteration0_RJMCMC scenario", {
  # Matching the exact scenario from iteration0_RJMCMC
  set.seed(123)
  T_seg4 <- cumsum(rexp(16, rate = 1))  # 16 random time points
  Bvec4 <- c(0.0, 2.0, 4.0, 6.0, 8.0, 10.0)  # 6 breakpoints (S+1 breakpoints, S=5 segments)
  S4 <- length(Bvec4) - 1  # S = 5
  
  cat("\n=== Test 3: Exact scenario from iteration0_RJMCMC ===\n")
  cat("T_seg4 size:", length(T_seg4), "\n")
  cat("Bvec4 size:", length(Bvec4), "\n")
  cat("S (number of segments):", S4, "\n")
  cat("T_seg4 range:", min(T_seg4), "to", max(T_seg4), "\n")
  cat("Bvec4 range:", min(Bvec4), "to", max(Bvec4), "\n")
  
  result4 <- cpp_findInterval(T_seg4, Bvec4)
  
  cat("Result4 length:", length(result4), "\n")
  cat("Min result4 value:", min(result4), "\n")
  cat("Max result4 value:", max(result4), "\n")
  cat("Unique values in result4:", unique(sort(result4)), "\n")
  
  # Check that values are in valid range [0, S] for make_table
  expect_true(all(result4 >= 0))
  expect_true(all(result4 <= S4))  # Should be <= S, not S+1
  expect_true(all(result4 < length(Bvec4)))  # Should be < breaks.size()
  
  if(any(result4 > S4)){
    cat("✗ PROBLEM FOUND: Result contains values > S!\n")
    cat("  Invalid values (for make_table):", result4[result4 > S4], "\n")
    cat("  S =", S4, ", but got values up to", max(result4), "\n")
  } else {
    cat("✓ All values are <= S, safe for make_table\n")
  }
  
  cat("✓ Test 3 PASSED\n")
})

test_that("cpp_findInterval - Test Case 4: Verify fix caps values at S (not S+1)", {
  # Test that cpp_findInterval now caps values at breaks.size() - 1 (S) instead of breaks.size() (S+1)
  # This test uses T_seg values > max(Bvec) to trigger the edge case
  T_seg5 <- c(seq(0.1, 4.9, length.out = 15), 6.0)  # Last element > max(Bvec)
  Bvec5 <- c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0)  # size 6 (S+1 breakpoints, S=5 segments)
  S5 <- length(Bvec5) - 1  # S = 5
  
  cat("\n=== Test 4: Verifying fix caps values at S (not S+1) ===\n")
  cat("T_seg5 size:", length(T_seg5), "\n")
  cat("Bvec5 size:", length(Bvec5), "\n")
  cat("S (number of segments):", S5, "\n")
  cat("T_seg5 range:", min(T_seg5), "to", max(T_seg5), "\n")
  cat("Bvec5 range:", min(Bvec5), "to", max(Bvec5), "\n")
  
  result5 <- cpp_findInterval(T_seg5, Bvec5)
  
  cat("Max result5:", max(result5), "\n")
  cat("Expected max (S):", S5, "\n")
  cat("Bvec5.size() (S+1):", length(Bvec5), "\n")
  
  # After the fix, max result should be <= S (i.e., <= breaks.size() - 1)
  # It should NOT equal breaks.size()
  expect_true(all(result5 >= 0))
  expect_true(all(result5 <= S5))  # Should be capped at S, not S+1
  expect_true(max(result5) < length(Bvec5))  # Should be strictly less than breaks.size()
  
  if(max(result5) == S5){
    cat("✓ Fix working: max result is capped at S =", S5, "\n")
    cat("  Values > max(Bvec) are correctly assigned to the last segment (index", S5, ")\n")
  } else if(max(result5) < S5){
    cat("✓ All values are within valid segment range [0,", S5, "]\n")
  } else {
    cat("✗ PROBLEM: max(result5) =", max(result5), "> S =", S5, "\n")
    stop("cpp_findInterval fix not working - still returns values > S!")
  }
  
  cat("✓ Test 4 PASSED - Fix verified\n")
})

