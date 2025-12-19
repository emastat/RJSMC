ts_data <- create_test_ts_data(n_obs = 150, num_logs = 5, seed = 1001)
parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)

# Define interval
start_point <- 3
end_point <- 10.0
t_star <- 2

# Filter observations within the interval
in_interval <- (ts_data$Tvec >= t_star) & (ts_data$Tvec < end_point)
T_seg <- ts_data$Tvec[in_interval]
Y_seg <- ts_data$Yvec[in_interval]

# Only proceed if we have observations

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
container_B <- vector("list", n_particle)
container_V <- vector("list", n_particle)
container_Z <- vector("list", n_particle)
container_Q <- vector("list", n_particle)
container_F <- vector("list", n_particle)
Svec <- integer(n_particle)
V_last_complete <- integer(n_particle)
#B_last <- numeric(n_particle)
B_last <- rep(1, n_particle)
weight_vec <- rep(1.0 / n_particle, n_particle)
particle_index_vec <- 0:(n_particle - 1)  # 0-indexed for C++

set.seed(2001)

results <- RJMCMC_SMC(
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
      )






##  TEST SMC TURCOTTE
ts_data <- create_test_ts_data(n_obs = 50, num_logs = 5, seed = 123)
parameters <- create_test_parameters(U = 3, W = 2, K = 2, num_logs = 5)

result <- SMC_turcotte_cpp(
      ts_data$Yvec, ts_data$Tvec,
      length_UI = 5.0,
      n_particle = 10,
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


# Define interval
start_point <- 0.0
end_point <- 10.0
sample_size <- 5

# Create mock breakpoint_list
# Each vector represents a particle's Bvec:
# - First element can be < start_point (carry-over from previous UI)
# - Middle elements are between start_point and end_point
# - Last element is > end_point (open segment boundary)
breakpoint_list <- list(
  c(0.5, 2.0, 5.0, 8.0, 12.0),    # Particle 1: Bvec[0]=-0.5, internal=[2,5,8], last=12
  c(1.0, 1.5, 4.5, 7.5, 11.0),   # Particle 2
  c(0.3, 2.5, 6.0, 9.0, 13.0),   # Particle 3
  c(0.8, 1.0, 3.0, 5.0, 7.0, 12.0), # Particle 4: more breakpoints
  c(0.2, 2.0, 4.0, 6.0, 8.0, 11.5)  # Particle 5
)

# Run the function
set.seed(123)
result <- breakpoints_sampling(
  start_point = start_point,
  end_point = end_point,
  breakpoint_list = breakpoint_list,
  sample_size = sample_size
)


english_words <- load_english_words_data()
  
  
n_particle <- 200

result <- SMC_turcotte_cpp(
  english_words$Yvec,
  english_words$Tvec,
  length_UI = 1.5,
  n_particle = n_particle,
  U = english_words$U,
  W = english_words$W,
  K = english_words$K,
  num_logs = english_words$num_logs,
  lambdamat = english_words$lambdamat,
  keyvec = english_words$keyvec,
  etavec = english_words$etavec,
  key0vec = english_words$key0vec,
  eta0vec = english_words$eta0vec,
  alphavec = english_words$alphavec,
  muvec = english_words$muvec,
  probvec_V = english_words$probvec_V,
  probvec_Z = english_words$probvec_Z,
  probvec_Q = english_words$probvec_Q,
  probvec_F = english_words$probvec_F,
  P0 = english_words$P0,
  minimum_n = english_words$minimum_n,
  Jss1 = 1/3,
  Js1s = 1/3,
  Smax = 150,
  n_ite = 20000,
  burn_in = 5000,
  thinning = 5
)


n_particle <- 300
  
result <- SMC_turcotte_cpp(
  real_data$Yvec[1:100],
  real_data$Tvec[1:100],
  length_UI = 2,
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
  n_ite = 2000,
  burn_in = 500,
  thinning = 5
)
