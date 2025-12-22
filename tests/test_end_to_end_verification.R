# End-to-End Test for RJSMC Algorithm Verification
# This test generates synthetic data with known true values for all state variables
# and breakpoints, then verifies that the SMC algorithm correctly approximates
# the true posterior distributions using visualization plots.

# Clear environment
rm(list=ls())

# Set seed for reproducibility
set.seed(12345)

# ============================================================================
# STEP 1: Define Model Parameters
# ============================================================================

# Number of states
K <- 2  # Number of Z states (length of non-empty segments)
U <- 4  # Number of V states (message composition)
W <- 2  # Number of Q states (rate of messages)

# Create lambda matrix for message composition (V states)
# Each row represents the probability distribution over LETTERS for each V state
lambdamat <- matrix(0, U, length(LETTERS))

# V state 1: "HELLO" pattern
lambdamat[1,] <- compute_lambdavec(c("H","E","L","L","O"), rep(0.85, 5)/5)

# V state 2: "WORLD" pattern  
lambdamat[2,] <- compute_lambdavec(c("W","O","R","L","D"), rep(0.85, 5)/5)

# V state 3: "TEST" pattern
lambdamat[3,] <- compute_lambdavec(c("T","E","S","T"), rep(0.85, 4)/4)

# V state 4: "DATA" pattern
lambdamat[4,] <- compute_lambdavec(c("D","A","T","A"), rep(0.85, 4)/4)

# Segment length parameters (Z states)
# Z=1: Long segments (mean 30 minutes, variance 5 minutes)
mu_L1 <- 30/60  # hours
var_L1 <- 5/60  # hours
keyvec <- c((mu_L1^2)/var_L1, (mu_L1^2)/var_L1 * 0.5)  # Z=2 is shorter
etavec <- c(mu_L1, mu_L1 * 0.5)

# Message rate parameters (Q states)
# Q=1: Low rate (60 messages/hour)
mu_m1 <- 60
var_m1 <- 10
# Q=2: High rate (300 messages/hour)
mu_m2 <- 300
var_m2 <- 30

alphavec <- c((mu_m1^2)/var_m1, (mu_m2^2)/var_m2)
muvec <- c(mu_m1, mu_m2)

# Empty segment length parameters (F states)
# F=1: Short empty periods (0.5 hours)
mu_Le1 <- 0.5
var_Le1 <- 0.1
# F=2: Long empty periods (2 hours)
mu_Le2 <- 2.0
var_Le2 <- 0.3

key0vec <- c((mu_Le1^2)/var_Le1, (mu_Le2^2)/var_Le2)
eta0vec <- c(mu_Le1, mu_Le2)

# Probability vectors for state transitions
probvec_V <- c(0.25, 0.25, 0.25, 0.25)  # Uniform for V states
probvec_Z <- c(0.6, 0.4)  # Prefer longer segments
probvec_Q <- c(0.7, 0.3)  # Prefer lower rates
probvec_F <- c(0.6, 0.4)  # Prefer shorter empty periods

P0 <- 0.3  # Probability of empty segment after non-empty

# ============================================================================
# STEP 2: Generate Synthetic Data with Known True Values
# ============================================================================

# Generate data using data_simulation function
# This creates coherent true values for Bvec, Vvec, Zvec, Qvec, Fvec
# and corresponding observations Tvec, Yvec
sim_data <- data_simulation(
  probvec_Z = probvec_Z,
  probvec_V = probvec_V,
  probvec_Q = probvec_Q,
  probvec_F = probvec_F,
  P0 = P0,
  alphavec = alphavec,
  muvec = muvec,
  keyvec = keyvec,
  etavec = etavec,
  key0vec = key0vec,
  eta0vec = eta0vec,
  lambdamat = lambdamat,
  K = K,
  U = U,
  W = W,
  Bmax = 25,  # Generate data up to time 25 hours
  seed = 12345,
  min_obs = 3
)

# Extract true values
true_Bvec <- sim_data$Bvec  # True breakpoints
true_Vvec <- sim_data$Vvec  # True V states (message composition)
true_Zvec <- sim_data$Zvec  # True Z states (segment length)
true_Qvec <- sim_data$Qvec  # True Q states (message rate)
true_Fvec <- sim_data$Fvec  # True F states (empty segment length)

# Extract observations
Tvec <- sim_data$Tvec  # Time stamps
Yvec <- sim_data$Yvec  # Message labels

# ============================================================================
# STEP 3: Prepare Data and Parameters for SMC
# ============================================================================

# Prepare time series data
ts_data <- list()
ts_data$Tvec <- Tvec
ts_data$Yvec <- Yvec

# Prepare parameters
parameters <- list()
parameters$U <- U
parameters$W <- W
parameters$K <- K
parameters$lambdamat <- lambdamat
parameters$keyvec <- keyvec
parameters$etavec <- etavec
parameters$key0vec <- key0vec
parameters$eta0vec <- eta0vec
parameters$alphavec <- alphavec
parameters$muvec <- muvec
parameters$probvec_V <- probvec_V
parameters$probvec_Z <- probvec_Z
parameters$probvec_Q <- probvec_Q
parameters$probvec_F <- probvec_F
parameters$P0 <- P0
parameters$minimum_n <- 0  # Allow segments with minimum observations

# Prepare settings
settings <- list()
settings$num_logs <- ncol(lambdamat)  # Number of possible messages (26 LETTERS)
settings$length_UI <- 2  # Length of each update interval (hours)
settings$n_particle <- 1000  # Number of particles
settings$Jss1 <- 1/3  # Probability to propose jump forward
settings$Js1s <- 1/3  # Probability to propose jump backward
settings$Smax <- 100  # Maximum number of segments per update interval
settings$n_ite <- 20000  # Number of RJMCMC iterations
settings$burn_in <- 5000  # Burn-in period
settings$thinning <- 5  # Thinning interval
settings$method <- "turcotte"  # Use Turcotte method

# Verify constraint: (n_ite - burn_in) / thinning > n_particle
if ((settings$n_ite - settings$burn_in) / settings$thinning <= settings$n_particle) {
  stop("Error: (n_ite - burn_in) / thinning must be > n_particle")
}

# ============================================================================
# STEP 4: Run SMC Algorithm
# ============================================================================

out_SMC <- RJSMC::SMC(
  ts_data = ts_data,
  parameters = parameters,
  settings = settings
)

# ============================================================================
# STEP 5: Visual Verification Using Plots
# ============================================================================

plot(out_SMC, 
     truth=list(B=true_Bvec, cl=true_Vvec),
     observations=list(Tvec=Tvec))
