# require("wasserstein_distance.R")
# require("charge_function_q.R")

##########
### 1D ###
##########

crit_1atatime = function(D, N, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != length(D)) stop("N is not the same as length of D")
  if(log_space == FALSE) {
    numPairs = N * (N - 1) / 2
    pairwise_PEs = rep(NA, numPairs)
    counter = 1
    qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p), D)
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_PEs[counter] = (qD[i] * qD[j] / sqrt((D[i] - D[j])^2))^k
        counter = counter + 1
      }
    }
    return((sum(pairwise_PEs))^(1/k))
  } else{
    if(N != length(D)) stop("N is not the same as length of D")
    numPairs = N * (N - 1) / 2
    pairwise_terms = rep(NA, numPairs)
    counter = 1
    logqD = sapply(FUN = function(x) log(q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                           f0, f1, type, var_margy0, var_margy1, p)), D)
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_terms[counter] = k * logqD[i] + k * logqD[j] - k * log(sqrt((D[i] - D[j])^2))
        counter = counter + 1
      }
    }
    return(exp((1 / k) * logSumExp(pairwise_terms)))
  }
}

##########
### 2D ###
##########

crit_1atatime_2d = function(D, N, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  # if(log_space == FALSE) {
  # numPairs = N * (N - 1) / 2
  # pairwise_PEs = rep(NA, numPairs)
  # counter = 1
  # qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  # for(i in 1:(N - 1)){
  #   for(j in (i + 1):N){
  #     pairwise_PEs[counter] = (qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2)))^k
  #     counter = counter + 1
  #   }
  # }
  # return((sum(pairwise_PEs))^(1/k))
  # } else{
  numPairs = N * (N - 1) / 2
  pairwise_terms = rep(NA, numPairs)
  counter = 1
  logqD = apply(D, 1, FUN = function(x) log(q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                                 f0, f1, type, var_margy0, var_margy1, p)))
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_terms[counter] = k * logqD[i] + k * logqD[j] - k * log(sqrt(sum((D[i, ] - D[j, ])^2)))
      counter = counter + 1
    }
  }
  return(exp((1 / k) * logSumExp(pairwise_terms)))
  # }
}

