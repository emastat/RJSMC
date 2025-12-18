# Tests for RJMCMC_SMC function

context("RJMCMC_SMC")

library(mclust)
library(fitdistrplus)

test_that("RJMCMC_SMC runs without error for a single interval", {
  # Test multiple scenarios with different parameters
  n_tests <- 20
  
  # Define test scenarios
  test_scenarios <- expand.grid(
    seed = 1001 + 0:(n_tests - 1),
    n_obs = c(30, 50, 100),
    t_star = c(0.0, 2.0, 5.0),
    start_point_base = c(0.0, 5.0, 10.0),
    interval_length = c(5.0, 10.0, 15.0)
  )
  
  # Limit to n_tests scenarios
  if (nrow(test_scenarios) > n_tests) {
    test_scenarios <- test_scenarios[sample(nrow(test_scenarios), n_tests), ]
  }
  
  for (test_case in 1:nrow(test_scenarios)) {
    scenario <- test_scenarios[test_case, ]
    
    # Create test data with different seeds and observation counts
    ts_data <- create_test_ts_data(
      n_obs = scenario$n_obs, 
      num_logs = 5, 
      seed = scenario$seed
    )
    parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)
    
    # Define interval with different start and end points
    start_point <- scenario$start_point_base
    end_point <- scenario$start_point_base + scenario$interval_length
    t_star <- scenario$t_star
    
    # Filter observations within the interval (matching SMC_turcotte pattern)
    # Use >= start_point and < end_point to match the C++ code
    in_interval <- (ts_data$Tvec >= start_point) & (ts_data$Tvec < end_point)
    T_seg <- ts_data$Tvec[in_interval]
    Y_seg <- ts_data$Yvec[in_interval]
    
    # Skip this test case if no observations
    if (length(T_seg) == 0) {
      next
    }
    
    # Ensure we have enough observations
    expect_true(length(T_seg) > 0, 
                info = paste("Test case", test_case, 
                             "- Must have observations in interval",
                             "- start_point =", start_point,
                             "- end_point =", end_point,
                             "- T_seg length =", length(T_seg)))
  
    # Set up RJMCMC parameters
    n_particle <- 10
    n_ite <- 1000
    burn_in <- 2
    thinning <- 1
    Smax <- 50
    minimum_n <- 1
    Jss1 <- 0.5
    Js1s <- 0.5
    empty_seg <- FALSE
    
    # Initialize output containers (modified by reference)
    # CRITICAL: These must be initialized exactly as in SMC_turcotte.cpp
    container_B <- vector("list", n_particle)
    container_V <- vector("list", n_particle)
    container_Z <- vector("list", n_particle)
    container_Q <- vector("list", n_particle)
    container_F <- vector("list", n_particle)
    
    # IntegerVector Svec(n_particle) in C++ creates vector of zeros
    Svec <- integer(n_particle)
    
    # V_last_complete and B_last should be initialized (from previous iteration in real usage)
    # In C++: IntegerVector V_last_complete = rep(1,n_particle)
    # In C++: NumericVector B_last = rep(0.0,n_particle)
    V_last_complete <- rep(1L, n_particle)  # Initialize to first state
    B_last <- rep(0.0, n_particle)          # Initialize to 0.0 (not start_point!)
    
    # NumericVector weight_vec = rep(1.0,n_particle) in C++
    weight_vec <- rep(1.0, n_particle)  # All weights = 1.0 initially (not normalized)
    
    # IntegerVector particle_index_vec - 0-indexed for C++
    # CRITICAL: Must be integer type for C++
    particle_index_vec <- as.integer(0:(n_particle - 1))  # 0-indexed for C++
    
    # Ensure particle_index_vec is not empty
    expect_true(length(particle_index_vec) > 0, 
                info = paste("Test case", test_case, "- particle_index_vec must not be empty"))
    
    # Call RJMCMC_SMC
    set.seed(2001 + test_case)  # Different seed for each test case
    expect_error(
      RJMCMC_SMC(
        T_seg = T_seg,
        Y_seg = Y_seg,
        U = parameters$U,
        K = parameters$K,
        W = parameters$W,
        start_point = start_point,
        end_point = end_point,
        t_star = t_star,
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
        minimum_n = minimum_n,
        Jss1 = Jss1,
        Js1s = Js1s,
        Smax = Smax,
        n_ite = n_ite,
        burn_in = burn_in,
        thinning = thinning,
        n_particle = n_particle,
        particle_index_vec = particle_index_vec,
        V_last_complete = V_last_complete,
        B_last = B_last,
        container_B = container_B,
        container_V = container_V,
        container_Z = container_Z,
        container_Q = container_Q,
        container_F = container_F,
        Svec = Svec,
        weight_vec = weight_vec,
        empty_seg = empty_seg
      ),
      NA,  # Expect no error
      info = paste("Test case", test_case,
                   "- RJMCMC_SMC should complete without error",
                   "- seed =", scenario$seed,
                   "- n_obs =", scenario$n_obs,
                   "- start_point =", start_point,
                   "- end_point =", end_point,
                   "- t_star =", t_star)
    )
    
    # Check that containers were modified
    expect_equal(length(container_B), n_particle,
                 info = paste("Test case", test_case, "- container_B length"))
    expect_equal(length(container_V), n_particle,
                 info = paste("Test case", test_case, "- container_V length"))
    expect_equal(length(container_Z), n_particle,
                 info = paste("Test case", test_case, "- container_Z length"))
    expect_equal(length(container_Q), n_particle,
                 info = paste("Test case", test_case, "- container_Q length"))
    expect_equal(length(container_F), n_particle,
                 info = paste("Test case", test_case, "- container_F length"))
    expect_equal(length(Svec), n_particle,
                 info = paste("Test case", test_case, "- Svec length"))
    expect_equal(length(V_last_complete), n_particle,
                 info = paste("Test case", test_case, "- V_last_complete length"))
    expect_equal(length(B_last), n_particle,
                 info = paste("Test case", test_case, "- B_last length"))
    expect_equal(length(weight_vec), n_particle,
                 info = paste("Test case", test_case, "- weight_vec length"))
    
    # Check that Svec contains valid values
    expect_true(all(Svec > 0),
                info = paste("Test case", test_case,
                             "- All Svec values should be > 0",
                             "- Svec =", paste(Svec, collapse = ", ")))
    expect_true(all(Svec <= Smax),
                info = paste("Test case", test_case,
                             "- All Svec values should be <= Smax",
                             "- Svec =", paste(Svec, collapse = ", "),
                             "- Smax =", Smax))
    
    # Check that containers are populated
    for (i in 1:n_particle) {
      expect_false(is.null(container_B[[i]]),
                   info = paste("Test case", test_case, "- container_B[[", i, "]] should not be NULL"))
      expect_false(is.null(container_V[[i]]),
                   info = paste("Test case", test_case, "- container_V[[", i, "]] should not be NULL"))
      expect_false(is.null(container_Z[[i]]),
                   info = paste("Test case", test_case, "- container_Z[[", i, "]] should not be NULL"))
      expect_false(is.null(container_Q[[i]]),
                   info = paste("Test case", test_case, "- container_Q[[", i, "]] should not be NULL"))
      expect_false(is.null(container_F[[i]]),
                   info = paste("Test case", test_case, "- container_F[[", i, "]] should not be NULL"))
      
      # Check dimensions match Svec[i]
      S_i <- Svec[i]
      expect_equal(length(container_B[[i]]), S_i + 1,
                   info = paste("Test case", test_case, "- container_B[[", i, "]] length should be", S_i + 1))
      expect_equal(length(container_V[[i]]), S_i,
                   info = paste("Test case", test_case, "- container_V[[", i, "]] length should be", S_i))
      expect_equal(length(container_Z[[i]]), S_i,
                   info = paste("Test case", test_case, "- container_Z[[", i, "]] length should be", S_i))
      expect_equal(length(container_Q[[i]]), S_i,
                   info = paste("Test case", test_case, "- container_Q[[", i, "]] length should be", S_i))
      expect_equal(length(container_F[[i]]), S_i,
                   info = paste("Test case", test_case, "- container_F[[", i, "]] length should be", S_i))
    }
  }
})

