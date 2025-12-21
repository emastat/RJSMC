library(RJSMC)
library(mclust)
library(fitdistrplus)

# Load real_data from package
data("real_data", package = "RJSMC")

set.seed(22)

# Prepare time series data
ts_data <- list()
ts_data$Yvec <- real_data$Yvec[1:1000]
ts_data$Tvec <- real_data$Tvec[1:1000]

# Prepare parameters from real_data
parameters <- list()
parameters$U <- real_data$U
parameters$W <- real_data$W
parameters$K <- real_data$K
parameters$lambdamat <- real_data$lambdamat
parameters$keyvec <- real_data$keyvec
parameters$etavec <- real_data$etavec
parameters$key0vec <- real_data$key0vec
parameters$eta0vec <- real_data$eta0vec
parameters$alphavec <- real_data$alphavec
parameters$muvec <- real_data$muvec
parameters$probvec_V <- real_data$probvec_V
parameters$probvec_Z <- real_data$probvec_Z
parameters$probvec_Q <- real_data$probvec_Q
parameters$probvec_F <- real_data$probvec_F
parameters$P0 <- real_data$P0
parameters$minimum_n <- real_data$minimum_n

# Prepare settings
settings <- list()
settings$num_logs <- real_data$num_logs
settings$length_UI <- 2
settings$n_particle <- 1500
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 50
settings$n_ite <- 20000
settings$burn_in <- 5000
settings$thinning <- 5
settings$method <- "turcotte"
settings$recycled_particles <- 4

# Run SMC end-to-end
cat("Starting SMC computation...\n")
cat("Data points:", length(ts_data$Yvec), "\n")
cat("Number of particles:", settings$n_particle, "\n")
cat("Update interval length:", settings$length_UI, "\n")
cat("Number of iterations:", settings$n_ite, "\n")
cat("\n")

out_SMC <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)

# Plot results
#plot(out_SMC, truth = NULL)

