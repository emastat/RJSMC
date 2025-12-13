#' @title compute_lambdavec
#' function for generating the lambda parameters, a matrix where each row is the probability distribution
#' for the log messages (here we use LETTERS) for each V state assumed
#'
#' @param name vector containing the LETTTERS we want to assign a specific probability
#' @param prob_name vector with assigned probabilities for the LETTERS specified in \code{name}
#' @return matrix with probability vector for each LETTER, in each V state
#' @export

compute_lambdavec <- function(name,prob_name){

  prob_vec <- rep(0,length(LETTERS))

  for(i in 1:length(name)){

    prob_vec[which(name[i]==LETTERS)] =prob_vec[which(name[i]==LETTERS)] + prob_name[i] }

  outside <- setdiff(LETTERS,name)

  for(j in 1:length(outside)){

    prob_vec[which(outside[j]==LETTERS)] = (1-sum(prob_name))/length(outside)

    }

  return(prob_vec) ;

  }
