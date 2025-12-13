library(RJSMC)
library(mclust)
library(fitdistrplus)

# load abb data and add the to the RJSMC package

abb_data <- readRDS("/Users/emanueleuio/Desktop/di21_final.RDS")
Tvec <- abb_data$time
Yvec <- as.numeric(abb_data$label)

set.seed(22)

ts_data <- list()
ts_data$Yvec <- Yvec[1:1000]
ts_data$Tvec <- Tvec[1:1000]
lambdamat <- readRDS("/Users/emanueleuio/Desktop/parameter_values/lambdamat_par.RDS")
keyvec <- readRDS("/Users/emanueleuio/Desktop/parameter_values/greenvec_par.RDS")
etavec <-  readRDS("/Users/emanueleuio/Desktop/parameter_values/etavec_par.RDS")
key0vec <- readRDS("/Users/emanueleuio/Desktop/parameter_values/greenvec_empty_par.RDS")
eta0vec <-  readRDS("/Users/emanueleuio/Desktop/parameter_values/etavec_empty_par.RDS")
alphavec <- readRDS("/Users/emanueleuio/Desktop/parameter_values/alphavec_par.RDS")
muvec <- readRDS("/Users/emanueleuio/Desktop/parameter_values/muvec_par.RDS")
probvec_Z <- readRDS("/Users/emanueleuio/Desktop/parameter_values/qvec_par.RDS")
probvec_V <- readRDS("/Users/emanueleuio/Desktop/parameter_values/deltavec_par.RDS")
probvec_F <- readRDS("/Users/emanueleuio/Desktop/parameter_values/epsivec_par.RDS")
probvec_Q <- c(0.5,0.5)
num_logs <- ncol(lambdamat)
P0 <- readRDS("/Users/emanueleuio/Desktop/parameter_values/P0_par.RDS")
minimum_n <- 0
probvec_Z <- readRDS("/Users/emanueleuio/Desktop/parameter_values/qvec_par.RDS")
probvec_V <- readRDS("/Users/emanueleuio/Desktop/parameter_values/deltavec_par.RDS")
probvec_F <- readRDS("/Users/emanueleuio/Desktop/parameter_values/epsivec_par.RDS")
probvec_Q <- c(0.5,0.5)

K <- length(alphavec)
U <- length(probvec_V)
W <- 2


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
parameters$minimum_n <- 0

settings <- list()
settings$num_logs <- num_logs
settings$length_UI <- 5
settings$n_particle <- 10000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 50
settings$n_ite <- 60000
settings$burn_in <- 10000
settings$thinning <- 5
settings$method <- "turcotte"
settings$recycled_particles = 4

#sink("/Users/emanueleuio/Desktop/RJMCMC-output.txt")


out_SMC <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)


plot(out_SMC, truth = NULL)

