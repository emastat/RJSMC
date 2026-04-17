# Official Simulation (LIGHT): SMC with different Update Interval lengths
#
# Same design as official_simulation.R (4 runs with length_UI = 0.5, 1, 2, 5)
# but with a lighter configuration for faster runs.
#
# Configuration (LIGHT):
#   - 1000 particles (vs 5000 in full)
#   - 10000 RJMCMC iterations (vs 50000 in full)
#   - 4000 burn-in iterations (vs 20000 in full)
#   - thinning = 5 (same as full)
#
# Report is written to: results_explanation/official_simulation_light.md
# (Different from official_simulation.md so full and light reports do not overwrite each other.)
#
# Full run: tests/official_simulation.R → official_simulation.md

# Clear environment
rm(list=ls())
set.seed(32453)

# Load required libraries
devtools::load_all()
library(mclust)
library(fitdistrplus)

# Load the "english_words" data
Bvec <- english_words$Bvec
Vvec <- english_words$Vvec
Zvec <- english_words$Zvec
Qvec <- english_words$Qvec
Fvec <- english_words$Fvec

Tvec <- english_words$Tvec
Yvec <- english_words$Yvec

# Prepare time series data
ts_data <- list()
ts_data$Yvec <- Yvec
ts_data$Tvec <- Tvec

# Prepare parameters (same for all runs)
parameters <- list()
parameters$U <- english_words$U
parameters$W <- english_words$W
parameters$K <- english_words$K
parameters$lambdamat <- english_words$lambdamat
parameters$keyvec <- english_words$keyvec
parameters$etavec <- english_words$etavec
parameters$key0vec <- english_words$key0vec
parameters$eta0vec <- english_words$eta0vec
parameters$alphavec <- english_words$alphavec
parameters$muvec <- english_words$muvec
parameters$probvec_V <- english_words$probvec_V
parameters$probvec_Z <- english_words$probvec_Z
parameters$probvec_Q <- english_words$probvec_Q
parameters$probvec_F <- english_words$probvec_F
parameters$P0 <- english_words$P0
parameters$minimum_n <- english_words$minimum_n

# Prepare base settings (LIGHT: fewer particles and iterations; length_UI will vary)
base_settings <- list()
base_settings$num_logs <- english_words$num_logs
base_settings$n_particle <- 1000
base_settings$Jss1 <- 1/3
base_settings$Js1s <- 1/3
base_settings$Smax <- 150
base_settings$n_ite <- 10000
base_settings$burn_in <- 4000
base_settings$thinning <- 5
base_settings$method <- "turcotte"


settings <- base_settings
settings$length_UI <- 2

#We need that (n_ite-burn_in)/thinning > n_particle

# SMC() returns raw C++ output; smc_post_processing() builds the RJSMC object for plot/breakpoint_cx
out_cpp <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)

disc <- seq(min(ts_data$Tvec), max(ts_data$Tvec), by = 0.25)  # or c(...), sort(unique(...)), etc.

out_SMC <- RJSMC::smc_post_processing(out_cpp,
                                      parameters = parameters,
                                      settings = settings,
                                      interval_length = 0.1,
                                      discretization_points = disc)

## plot with observations
plot(out_SMC, 
     truth=list(B=english_words$Bvec, cl=english_words$Vvec),
     observations=list(Tvec=Tvec),time_to_date=TRUE)


plot(out_SMC, 
     truth=list(B=english_words$Bvec, cl=english_words$Vvec),
     observations=list(Tvec=Tvec),time_to_date=TRUE, pl="B")

plot(out_SMC, 
     truth=list(B=english_words$Bvec, cl=english_words$Zvec),
     observations=list(Tvec=Tvec),time_to_date=TRUE, pl="Z")     

#  Test breakpoint_cx
test = breakpoint_cx(
     out_SMC,
     true_Bvec=english_words$Bvec,
     Tvec=Tvec,
     U=english_words$U,
     storage_weight=out_SMC@storage_weight,
     UI_bounds=out_SMC@UI_bounds,
     Vvec=english_words$Vvec)

