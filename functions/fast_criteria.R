##########
### 1D ###
##########

crit_fast = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / sqrt((D[i] - D[j])^2)
      counter = counter + 1
    }
  }
  return(max(pairwise_PEs))
}

##########
### 2D ###
##########

crit_fast_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  # if(log_space == FALSE){
  #   numPairs = N * (N - 1) / 2
  #   pairwise_PEs = rep(NA, numPairs)
  #   counter = 1
  #   qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  #   for(i in 1:(N - 1)){
  #     for(j in (i + 1):N){
  #       pairwise_PEs[counter] = qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2))
  #       counter = counter + 1
  #     }
  #   }
  #   return(max(pairwise_PEs))
  # } else{
  numPairs = N * (N - 1) / 2
  pairwise_terms = rep(NA, numPairs)
  counter = 1
  logqD = apply(D, 1, FUN = function(x) log(q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                                 f0, f1, type, var_margy0, var_margy1, p)))
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_terms[counter] = logqD[i] + logqD[j] - log(sqrt(sum((D[i, ] - D[j, ])^2)))
      counter = counter + 1
    }
  }
  return(max(exp(pairwise_terms)))
  # }
}