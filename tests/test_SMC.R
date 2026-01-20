#  test SMC.R

#load the "english_words" data
#data("english_words")
rm(list=ls())
set.seed(32453)    ##Works
#set.seed(45635)     ##Gets an error

devtools::load_all()

library(mclust)
library(fitdistrplus)

Bvec <- english_words$Bvec
Vvec <- english_words$Vvec
Zvec <- english_words$Zvec
Qvec <- english_words$Qvec
Fvec <- english_words$Fvec

Tvec <- english_words$Tvec
Yvec <- english_words$Yvec

ts_data <- list()
ts_data$Yvec <- Yvec
ts_data$Tvec <- Tvec

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

settings <- list()
settings$num_logs <- english_words$num_logs
settings$length_UI <- 2
settings$n_particle <- 5000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 150
settings$n_ite <- 50000
settings$burn_in <- 20000
settings$thinning <- 5
settings$method <- "turcotte"

#We need that (n_ite-burn_in)/thinning > n_particle

out_SMC <- RJSMC::SMC(ts_data = ts_data,
                    parameters = parameters,
                    settings = settings)

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

