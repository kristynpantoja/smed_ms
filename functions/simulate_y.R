simulateY = function(D, N, mean_beta, var_mean, var_e, numSims, type = NULL, seed = NULL){
  if(is.null(seed)) set.seed(123)
  X = constructDesignX(D, N, type)
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    beta = t(rmvnorm(n = 1, mean = mean_beta, sigma = var_mean))
    for(i in 1:N){
      Y[i, j] = rnorm(n = 1, mean = X[i, ] %*% beta, sd = sqrt(var_e))
    }
  }
  return(Y)
}