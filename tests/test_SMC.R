#  test SMC.R

#load the "english_words" data
#data("english_words")
rm(list=ls())
#set.seed(32453)    ##Works
set.seed(45635)     ##Gets an error
library(RJSMC)
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
settings$length_UI <- 0.5
settings$n_particle <- 5000
settings$Jss1 <- 1/3
settings$Js1s <- 1/3
settings$Smax <- 150
settings$n_ite <- 40000
settings$burn_in <- 6000
settings$thinning <- 5
settings$method <- "turcotte"

#We need that (n_ite-burn_in)/thinning > n_particle

out_SMC <- RJSMC::SMC(ts_data = ts_data,
                    parameters = parameters,
                    settings = settings)


plot(out_SMC,truth=list(B=english_words$Bvec,cl=english_words$Vvec))

par1 = english_words$alphavec[1]
par2 = (english_words$alphavec[1]/english_words$muvec)/(english_words$alphavec[1]/english_words$muvec+ 15.1652-14.495)

dnbinom(-1,par1,par2,log = TRUE )
# # bugfix get_results
#
# UI_bounds= readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/UI_bounds.rds")
# n_particle = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/n_particle.rds")
# storage_B = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/out_SMC_cpp$storage_B.rds")
# storage_V = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/out_SMC_cpp$storage_V.rds")
# storage_Z = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/out_SMC_cpp$storage_Z.rds")
# storage_Q = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/out_SMC_cpp$storage_Q.rds")
# storage_F = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/out_SMC_cpp$storage_F.rds")
# U = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/U.rds")
# W = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/W.rds")
# K = readRDS(file ="/Users/emanueleuio/Desktop/inspect_results/K.rds")
#
# j = 1
# results_UI <- get_results(UI_bounds[j],
#                           UI_bounds[j+1],
#                           n_particle,
#                           0.01,
#                           storage_B[[j]],
#                           storage_V[[j]],
#                           storage_Z[[j]],
#                           storage_Q[[j]],
#                           storage_F[[j]],
#                           U,
#                           W,
#                           K,
#                           2)
#
#
#
#
