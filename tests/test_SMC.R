#  test SMC.R

#load the "english_words" data
#data("english_words")
library(RJSMC)

Bvec <- english_words$Bvec
Vvec <- english_words$Vvec
Zvec <- english_words$Zvec
Qvec <- english_words$Qvec
Fvec <- english_words$Fvec

Tvec <- english_words$Tvec
Yvec <- english_words$Yvec

ts_data <- list()
#ts_data$Yvec <- Yvec[Tvec>Bvec[1] & Tvec<Bvec[20]]
#ts_data$Tvec <- Tvec[Tvec>Bvec[1] & Tvec<Bvec[20]]

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
settings$length_UI <- 4
settings$n_particle <- 1000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 150
settings$n_ite <- 20000
settings$burn_in <- 10000
settings$thinning <- 10
settings$method <- "waste_free"


out_SMC <- RJSMC::SMC(ts_data = ts_data,
                      parameters = parameters,
                      settings = settings)


plot(out_SMC,truth=list(B=english_words$Bvec,cl=english_words$Vvec))


