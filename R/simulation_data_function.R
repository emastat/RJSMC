#' @title function for generating simulation data
#'  This function, given a set of hyper-parameters, generates a complete dataset. The generating
#' model is the model of paper 1 but assuming 3 states (Z,Q,V) for length, rate and composition
#' of the segments, respectively.
#'
#' @param probvec_Z Probability vector for state Z
#' @param probvec_V Probability vector for state V
#' @param probvec_Q Probability vector for state Q
#' @param probvec_F Probability vector for state F
#' @param P0        Probability to observe an empty segment after a non-empty one
#' @param alphavec  Shape parameters for the rate of the non-empty segments
#' @param muvec     Mean parameters for the rate of the non-empty segments
#' @param  keyvec Shape parameters for the length of the non-empty segments
#' @param  etavec   Mean parameters for the length of the non-empty segments
#' @param key0vec Shape parameters for the length of the empty segments
#' @param eta0vec   Mean parameters for the length of the empty segments
#' @param lambdamat      Probability vectors for the message composition of the V states
#' @param K              Number of levels for state Z
#' @param U              Number of levels for state V
#' @param W              Number of levels for state Q
#' @param Bmax           Final breakpoint
#' @param seed           Seed for reproducibility
#' @param  min_obs       Minimum number of observations allowed in a segment, default is 3
#' @return list with following objects
#' \describe{
#' \item{Bvec}{ Vector with generated changepoints}
#' \item{Vvec}{Vector with generated V's - ruling the message distribution in the non-empty segments}
#' \item{Zvec}{Vector with generated Z's - ruling the length of  non-empty segments}
#' \item{Qvec}{Vector with generated Q's - ruling the rate of messages in non-empty segments}
#' \item{Fvec}{Vector with generated F's - ruling the length of the empty segments}
#' \item{Mvec}{Vector with generated M's - rates of messages in non-empty segments}
#' \item{Nvec}{Vector with generated N's - counts of observations in non-empty segments}
#' \item{Tvec}{Vector with genetated time stamps}
#' \item{Yvec}{Vector with genetated messages}
#' \item{probvec_Z}{vector with probabilities for the levels of state Z}
#' \item{probvec_V}{vector with probabilities for the levels of state V}
#' \item{probvec_F}{vector with probabilities for the levels of state F}
#' \item{probvec_Q}{vector with probabilities for the levels of state Q}
#' \item{P0}{Probability to observe an empty segment after a non-empty one}
#' \item{alphavec}{Shape parameters for the rate of the non-empty segments}
#' \item{muvec}{Mean parameters for the rate of the non-empty segments}
#' \item{keyvec}{Shape parameters for the length of the non-empty segments}
#' \item{key0vec}{Shape parameters for the lenght of the empty segments }
#' \item{etavec}{Mean parameters for the length of the non-empty segments}
#' \item{eta0vec}{ Mean parameters for the lenght of the empty segments}
#' \item{lambdamat}{Probability vectors for the message composition of the V states}
#' \item{K}{Number of levels for state Z}
#' \item{U}{Number of levels for state V}
#' \item{W}{Number of levels for state Q}
#' \item{num_logs}{number of different possible messages}
#' \item{minimum_n}{minimum number of messages allowed in a segment}
#'}
#' @export
data_simulation <- function(probvec_Z, #prob vector for state Z
                            probvec_V, #prob vector for state V
                            probvec_Q, #prob vector for state Q
                            probvec_F, #prob vector for the occurrence of state empty state
                            P0,        # parameter P0
                            alphavec, #shape parameter for intensity M
                            muvec,       # mean parameter for the intensity M
                            keyvec,     #shape parameter for  L
                            etavec,        #mean parameter for L
                            key0vec, # shape parameter for L in empty segments
                            eta0vec,   # mean parameter for L in empty segments
                            lambdamat,   #prob vector for the logs in each state V
                            K, # number of states Z
                            U, # number of states V
                            W, # number of states Q
                            Bmax, #final breakpoint
                            seed,  # seed
                          min_obs =3 # minimum number of obs in a segment
                            ){

  #1) set the seed
  set.seed(seed)

  #2) create a WHILE loop that generate at each iteration, a new interval and all the related quantities, until "Bmax" is reached

  Tj <- 0  # the first breakpoint
  V_old <- 1   # the initial V state. It's set to 1 so that the first segment can belong to any V state
  Bvec <- Tj    # initialize the breakpoint vector with Tj
  Zvec <- NULL  # vector with states ruling the length of the non-empty segments
  Qvec <- NULL  # vector with states ruling the rate of the non-empty segments
  Vvec <- NULL  # vector states ruling the message composition of the segments
  Mvec <- NULL  # vector hosting the rate in each segment
  Nvec <- NULL  # vector hosting the count in each segment
  Fvec <- NULL  # vector with state ruling the length of the empty segments
  time_stamps <- NULL # vector hosting the time-stamps observations
  log_events <- NULL  # vector hosting the messages observations
  V_new <- 0 # the new value for the V state of the segment

  while(Tj<Bmax){


    #generate the state V for interval L_j

    if(V_old==0){

      V_new <-  sample(x=1:U,size=1,prob=probvec_V )

    }else{

      V_new <-  sample(x=0:U,size=1,prob=c(P0, (1-P0)*probvec_V )  )

    }


    # update Vvec
    Vvec <- c(Vvec,V_new)


    if(V_new == 0){

      # generate length of the empty segment
      F_value <- sample(x=1:2,size = 1,prob = probvec_F)

      Fvec <- c(Fvec,F_value)

      L <- rgamma(n=1,shape=key0vec[F_value],rate=key0vec[F_value]/eta0vec[F_value])

      Tj <- Tj + L


      Bvec <- c( Bvec, Tj )

      # set to 0 the other states
      Zvec <- c(Zvec,0)
      Qvec <- c(Qvec,0)
      Mvec <- c(Mvec,0)
      Nvec <- c(Nvec,0)

    }else{

      #generate a non-empty segment


      # dample A and the the length
      Z <- sample(x=1:K,size=1,prob=probvec_Z)

      Zvec <- c(Zvec,Z)

      L <- rgamma(n=1,shape=keyvec[Z],rate=keyvec[Z]/etavec[Z])

      if(L<1/7200){L <- 1/7200}


     # sample Q and the the rate
      Q <- sample(x=1:W,size=1, prob = probvec_Q)

      Qvec <- c(Qvec,Q)

      m <- rgamma(n=1,shape=alphavec[Q],rate=alphavec[Q]/muvec[Q])

      Mvec <- c(Mvec,m)


      # set F to 0 because it's a non-empty segment
      Fvec <- c(Fvec,0)


      # generate the count of observations within the segment
      x <- rpois(n=1,lambda = L*m)

      x <- x + min_obs

      Nvec <- c(Nvec,x)



      # generate the position of the observations within the segment
      y <- sort(runif( n=x,min=Tj,max=Tj + L ))

      time_stamps <- c(time_stamps, y )

      e <- sample(x=1:ncol(lambdamat),size=x,prob =lambdamat[V_new,],replace = T )

      log_events <- c( log_events, e )

      Tj <- Tj + L


      Bvec <- c( Bvec, Tj )
    }

    V_old <- V_new

  }

  #remove the last breakpoint


  #STORE THE RESULTS

  output <- list()
  output$Bvec <- Bvec
  output$Vvec <- Vvec
  output$Zvec <- Zvec
  output$Qvec <- Qvec
  output$Fvec <- Fvec
  output$Mvec <- Mvec
  output$Nvec <- Nvec
  output$Tvec <- time_stamps
  output$Yvec <- log_events

  output$probvec_Z <- probvec_Z
  output$probvec_V <- probvec_V
  output$probvec_Q <- probvec_Q
  output$probvec_F <- probvec_F
  output$P0 <- P0
  output$alphavec <- alphavec
  output$muvec <- muvec
  output$keyvec <- keyvec
  output$key0vec <- key0vec
  output$etavec <- etavec
  output$eta0vec <- eta0vec
  output$lambdamat <- lambdamat

  output$K <- K
  output$U <- U
  output$W <- W
  output$minimum_n <- min_obs
  output$num_logs <- ncol(output$lambdamat)

  return(output)


}

