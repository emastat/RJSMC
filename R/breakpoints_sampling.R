#' @title  Function for generating a sample of Breakpoints inside an observational interval
#' Function for generating breakpoints from a kernel density inside an interval
#'
#' @param  start_point start point for the segment (t_star)
#' @param  end_point end point of the observational interval
#' @param breakpoint_list list with breakpoints vectors sampled during the RJSMCMC
#' @param sample_size number of samples to generate
#' @return list called output with following objects
#' \describe{
#' \item{logdens}{vector with log density value for  the generate breakpoint vector, for each sample }
#' \item{ S.samp}{list of vectors with the generated breakpoints in each sample }
#' @export
breakpoints_sampling = function(
    start_point,
    end_point,
    breakpoint_list,
    sample_size
)
{

  # unlist breakpoints in sorted vector and remove breakpoint smaller than start_point ( i.e t_star)

  copy_breakpoint_list = lapply(breakpoint_list,identity)

  breakpoints_flat = do.call(c,copy_breakpoint_list)
  sorted_breakpoints = sort(breakpoints_flat)
  jittered_breakpoints = jitter(sorted_breakpoints,amount = 0.01)
  jittered_breakpoints_inside = jittered_breakpoints[(jittered_breakpoints>start_point) &
                                            (jittered_breakpoints<end_point)]

  #and remove start_point
  Sv = jittered_breakpoints_inside - start_point

  if(length(Sv)==1){Sv = c(Sv,Sv)}


  if(any(Sv <=0)){

    stop("Breakpoint_sampling.R Breakpoints with value <=0 have been found. Breakpoint must be positive")
  }

  # length of the observational interval
  delta = end_point - start_point

  #  apply logit transformation. Allow breakpoint to be in (-inf, +inf)
  Sv.tr = log((Sv/delta)/(1-Sv/delta))

  if(any(is.nan(Sv.tr))){stop("breakpoint_sampling.R --> check line 45")}

  mclust_error = TRUE

  while(mclust_error){

    result <- try({
    fit =Mclust(Sv.tr)
                })

    if (class(result) == "try-error"){

      stop("breakpoint_sampling.R --> check line 57")

    }else{

      # fit MVN mixture to the breakpoints
      fit =Mclust(Sv.tr, verbose = FALSE)
      mclust_error = FALSE

    }
  }

  # simulate number of breakpoints from a negative binomial, based on the empirical counts
  # of breakpoints in the particles

  if(length(copy_breakpoint_list)<2){

    empirical_breakpoint_dist = rep(length(copy_breakpoint_list[1]),2)

  }else{empirical_breakpoint_dist = sapply(breakpoint_list, length)}


  suppressWarnings(fit.N <- fitdist(empirical_breakpoint_dist, "nbinom"))

  # sample number of breakpoint for each sample

  N.samp = rnbinom(sample_size, size=fit.N$estimate[1], mu=fit.N$estimate[2])

  #Simulate breakpoint positions for each the N.samp breakpoints simulated, forall the "sample_size" samples & compute the overall density


  S.samp  = vector(mode='list', length=length(N.samp))

  logdens = rep(NA,sample_size)

  for(i in 1:sample_size){

    logdens[i] = dnbinom(N.samp[i],size=fit.N$estimate[1],mu=fit.N$estimate[2],
                         log = TRUE)

    if(N.samp[i]>0){

      #Simulate  breakpoints on transformed scale

      #foo = sort(sim(fit$modelName,fit$parameters,n=N.samp[i])[,2])

      foo <- tryCatch({
      sort(sim(fit$modelName,fit$parameters,n=N.samp[i])[,2])
      },
      error = function(e) {
        # Code to handle the error
        stop("breakpoints_sampling.R-->: ", e$message)
      })


      #Transform back
      foo_exp = delta*exp(foo)/(1+exp(foo))

      #Calculate Jacobian on log-scale
      logJac = log(foo_exp)+log(delta-foo_exp)-log(delta)

      # update the log density
      logdens[i] = logdens[i] +
        sum(dens(foo,fit$modelName,fit$parameters,log=TRUE)+logJac)

      # numerical stability
      if (logdens[i]<-100) {

        logdens[i] = -100

      }

      S.samp[[i]] = unique(c(0, foo_exp, end_point-start_point))

      # START CODE TO DELETE
      if (any(is.nan(S.samp[[i]]))) {
        stop("breakpoints_sampling.R --> NaN value detected in the vector.")
      }
      #END CODE TO DELETE

    }else{

      foo = c(0, end_point-start_point)
      S.samp[[i]] = foo



    }
  }

  # add to each sampled vector start_point

  S.samp = lapply(S.samp, function(x) x+ start_point)


  # return list with sampled breakpoints and logdensity value

  output = list()

  output$log_density = logdens
  output$S.samp = S.samp

  return(output)


}
