# require("construct_design_matrix.R")

simulateY = function(D, N, true_beta, var_e, numSims, type = NULL, seed = NULL){
  # note that type is based on what type of model for true_beta
  if(is.null(seed)) set.seed(123)
  X = constructDesignX(D, N, type)
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    Y[ , j] = rmvnorm(1, X %*% true_beta, var_e * diag(N))
  }
  return(Y)
}



# simulateY_old = function(D, N, mean_beta, var_mean, var_e, numSims, type = NULL, seed = NULL){
#   if(is.null(seed)) set.seed(123)
#   X = constructDesignX(D, N, type)
#   Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
#   for(j in 1:numSims){
#     beta = t(rmvnorm(n = 1, mean = mean_beta, sigma = var_mean))
#     for(i in 1:N){
#       Y[i, j] = rnorm(n = 1, mean = X[i, ] %*% beta, sd = sqrt(var_e))
#     }
#   }
#   return(Y)
# }

