devtools::load_all()
library(mclust)
library(fitdistrplus)

# Load real_data from package
data("real_data", package = "RJSMC")

set.seed(22)

### Update interval 2.5 hours ###

# Prepare time series data
ts_data <- list()
ts_data$Yvec <- real_data$Yvec[1100:2200]
ts_data$Tvec <- real_data$Tvec[1100:2200]

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
settings$length_UI <- 2.5
settings$n_particle <- 1000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 50
settings$n_ite <- 10000
settings$burn_in <- 4000
settings$thinning <- 5
settings$method <- "turcotte"

out_SMC_2_5 <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)

# Plot results


png("real_data_state_reconstruction_2_5.png", width = 1600, height = 900, res = 200)

plot(
    out_SMC_2_5,
    observations=list(Tvec=ts_data$Tvec),
    time_to_date=FALSE,
    title= "State reconstruction - Update interval 2.5 hours")

dev.off()


########################################################

### Update interval 0.5 hours ###

# Prepare time series data
ts_data <- list()
ts_data$Yvec <- real_data$Yvec[1100:2200]
ts_data$Tvec <- real_data$Tvec[1100:2200]

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
settings$length_UI <- 0.5
settings$n_particle <- 1000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 50
settings$n_ite <- 10000
settings$burn_in <- 4000
settings$thinning <- 5
settings$method <- "turcotte"

out_SMC_0_5 <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)

# Plot results


png("real_data_state_reconstruction_0_5.png", width = 1600, height = 900, res = 200)

plot(
    out_SMC_0_5,
    observations=list(Tvec=ts_data$Tvec),
    time_to_date=FALSE,
    title= "State reconstruction - Update interval 0.5 hours")

dev.off()



########################################################

## COMP


cmp <- compare_smc_runs_js(
  obj_a = out_SMC_2_5,
  obj_b = out_SMC_0_5,
  Tvec = ts_data$Tvec,
  U = real_data$U,
  plot = TRUE,
  plot_path = "js_comparison_real_data.png",  # or any path you want
  plot_width = 7,
  plot_height = 4.5
)
