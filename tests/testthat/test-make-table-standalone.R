# Standalone test for make_table to identify out-of-bounds access
# This tests the function with interval values that could cause issues

context("make_table - Standalone Debug Test")

test_that("make_table with interval size 16 and S=15", {
  # Test Case: interval has 16 elements, S=15 (Bvec has 16 breakpoints)
  # This matches the scenario where we get "index 16 >= vector size 16"
  set.seed(123)
  interval <- sample(0:15, size = 16, replace = TRUE)  # 16 elements, values 0-15
  S <- 15  # 15 segments (Bvec would have 16 breakpoints)
  
  cat("\n=== Test 1: interval size 16, S=15 ===\n")
  cat("interval size:", length(interval), "\n")
  cat("S (number of segments):", S, "\n")
  cat("interval values:", sort(unique(interval)), "\n")
  cat("Max interval value:", max(interval), "\n")
  
  # make_table creates freque(K+1) where K=S, so freque has size S+1 = 16
  # It loops for(int i=0; i<K+1; i++) which is i=0 to i=15
  # Valid indices for freque are 0-15
  
  tryCatch({
    result <- make_table(interval, S, FALSE)
    
    cat("Result size:", length(result), "\n")
    cat("Expected size (S):", S, "\n")
    cat("Result values:", result, "\n")
    
    # After erasing index 0, result should have S elements
    expect_equal(length(result), S)
    expect_true(all(result >= 0))
    
    cat("✓ Test 1 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 1 FAILED:", e$message, "\n")
    stop(e)
  })
})

test_that("make_table with interval containing value = S", {
  # Test Case: interval contains value S (which should be valid)
  # make_table loops i=0 to i=S, so value S should be counted
  interval2 <- c(rep(0:14, each = 2), 15)  # Values 0-15, with 15 appearing once
  S2 <- 15
  
  cat("\n=== Test 2: interval contains value = S ===\n")
  cat("interval2 size:", length(interval2), "\n")
  cat("S2:", S2, "\n")
  cat("Max interval2 value:", max(interval2), "\n")
  cat("Unique values:", sort(unique(interval2)), "\n")
  
  tryCatch({
    result2 <- make_table(interval2, S2, FALSE)
    
    cat("Result2 size:", length(result2), "\n")
    cat("Result2 values:", result2, "\n")
    
    expect_equal(length(result2), S2)
    cat("✓ Test 2 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 2 FAILED:", e$message, "\n")
    stop(e)
  })
})

test_that("make_table with interval containing value > S (should be caught)", {
  # Test Case: interval contains value S+1 (which would be out of bounds)
  # This should not happen after our cpp_findInterval fix, but let's test
  interval3 <- c(rep(0:14, each = 2), 16)  # Values 0-14 and 16 (S+1)
  S3 <- 15
  
  cat("\n=== Test 3: interval contains value > S (S+1) ===\n")
  cat("interval3 size:", length(interval3), "\n")
  cat("S3:", S3, "\n")
  cat("Max interval3 value:", max(interval3), "\n")
  cat("Values > S3:", interval3[interval3 > S3], "\n")
  
  # make_table loops i=0 to i=S, so it never checks for value S+1
  # This means value S+1 won't be counted, but it shouldn't cause out-of-bounds
  # UNLESS make_table uses interval values as indices somewhere
  
  tryCatch({
    result3 <- make_table(interval3, S3, FALSE)
    
    cat("Result3 size:", length(result3), "\n")
    cat("Result3 values:", result3, "\n")
    
    # Value 16 won't be counted (loop only goes to 15)
    # But this shouldn't cause an error, just missing count
    expect_equal(length(result3), S3)
    cat("✓ Test 3 PASSED (value > S not counted, but no error)\n")
  }, error = function(e) {
    cat("✗ Test 3 FAILED:", e$message, "\n")
    cat("  This might indicate make_table uses interval values as indices!\n")
    stop(e)
  })
})

test_that("make_table exact scenario from iteration0_RJMCMC", {
  # Test Case: Exact scenario matching the error case
  # T_seg has 16 elements, Bvec has 16 breakpoints (S=15)
  set.seed(123)
  
  # Simulate what cpp_findInterval would return
  # After our fix, it should return values 0-15 (capped at S)
  interval4 <- sample(0:15, size = 16, replace = TRUE)
  S4 <- 15
  
  cat("\n=== Test 4: Exact scenario from iteration0_RJMCMC ===\n")
  cat("interval4 size:", length(interval4), "\n")
  cat("S4:", S4, "\n")
  cat("interval4 values:", sort(unique(interval4)), "\n")
  cat("Max interval4 value:", max(interval4), "\n")
  
  tryCatch({
    result4 <- make_table(interval4, S4, FALSE)
    
    cat("Result4 size:", length(result4), "\n")
    cat("Expected size (S4):", S4, "\n")
    
    expect_equal(length(result4), S4)
    expect_true(all(result4 >= 0))
    
    cat("✓ Test 4 PASSED\n")
  }, error = function(e) {
    cat("✗ Test 4 FAILED:", e$message, "\n")
    cat("  This is the exact scenario causing the out-of-bounds warning!\n")
    stop(e)
  })
})

