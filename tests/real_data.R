# load abb data and add the to the RJSMC package

abb_data <- readRDS("/Users/emanueleuio/Desktop/di21_final.RDS")
Tvec <- abb_data$time
Yvec <- as.numeric(abb_data$label)


ts_data <- list()
ts_data$Yvec <- Yvec[1:3000]
ts_data$Tvec <- Tvec[1:3000]
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
settings$length_UI <- 4
settings$n_particle <- 100
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 50
settings$n_ite <- 6000
settings$burn_in <- 2000
settings$thinning <- 10
settings$method <- "waste_free"
settings$recycled_particles = 4

#sink("/Users/emanueleuio/Desktop/RJMCMC-output.txt")


out_SMC <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)


# get number of Updating windows (i.e. intervals)
out_SMC@n_UI

# get marginal posterior for "V" at each discretization point

class(out_SMC@posteriors_container_V)
length(out_SMC@points_container)
table(out_SMC@UI_index_vector)


output =plot(out_SMC,
    truth = NULL)
output







