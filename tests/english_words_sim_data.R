library(extraDistr)

options(digits=9)


# number of Z and V states
K <- 2
U <- 5

lambdamat <- matrix(0,U,length(LETTERS))
lambdamat[1,] <-  compute_lambdavec(c("M","O","T","H","E","R"),rep(0.9,6)/6)
lambdamat[2,] <-  compute_lambdavec(c("B","E","A","C","H"),rep(0.9,5)/5)
lambdamat[3,] <-  compute_lambdavec(c("B","E","A","R"),rep(0.9,4)/4)
lambdamat[4,] <-  compute_lambdavec(c("F","A","T","H","E","R"),rep(0.9,6)/6)
lambdamat[5,] <-  compute_lambdavec(c("S","P","I","D","E","R"),rep(0.9,6)/6)


#20 minutes length (var 3 sec min)
mu_L1 <- 20/60 ; var_L1 <- 3/60

# 60 messages/hours : 60*(20/60) <- 20 each 20 minutes
mu_m1 <- 60; var_m1 <- 5

######

#5 minute length (var 40 sec)
mu_L2 <- 5/60 ; var_L2 <- 40/3600

# 2400  messages/hours : 2400/60 <- 40 each  minutes
mu_m2 <- 2400; var_m2 <- 20

#####
keyvec <- c((mu_L1^2)/var_L1,(mu_L2^2)/var_L2)

etavec <- c(mu_L1,mu_L2)

alphavec <- c((mu_m1^2)/var_m1,(mu_m2^2)/var_m2)

muvec <- c(mu_m1,mu_m2)

probvec_V <- c(0.13,0.3,0.15,0.2,0.22)

probvec_Z <- c(0.5,0.5)

probvec_Q <- c(0.7,0.3)

P0 <- 0.5

probvec_F <- c(0.4,0.6)

#empty vec 1 --> 1 hour length (var 1 min)
mu_Le1 <- 1 ; var_Le1 <- 6/60

##empty vec 2 --> 3 hours length (var 1 min)
mu_Le2 <- 3 ; var_Le2 <- 15/60

key0vec <- c((mu_Le1^2)/var_Le1,(mu_Le2^2)/var_Le2)

eta0vec <-c(mu_Le1,mu_Le2)

#<- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <-

english_words <- data_simulation(probvec_Z = probvec_Z,
                           probvec_V = probvec_V,
                           probvec_Q = probvec_Q,
                           probvec_F = probvec_F,
                           P0=P0,
                           alphavec=alphavec,
                           muvec=muvec,
                           keyvec = keyvec,
                           etavec= etavec,
                           key0vec = key0vec,
                           eta0vec = eta0vec,
                           lambdamat = lambdamat,
                           K = 2,
                           U = 5,
                           W = 2,
                           Bmax=50,
                           seed=52,
                           min_obs = 3)

# check if the generated datas  et is satisfactory


time_stamps <- english_words$Tvec
log_events <- english_words$Yvec
Bvec <- english_words$Bvec
Vvec <- english_words$Vvec
Zvec <- english_words$Zvec
Qvec <- english_words$Qvec
Fvec <- english_words$Fvec
Nvec <- english_words$Nvec

#add dataset to the RJSMC package
#usethis::use_data(english_words,RJSMC,overwrite = T)
