# require("construct_design_matrix.R")

simulateY = function(D, N, true_beta, var_e, numSims, type = NULL, seed = NULL){
  # note that type is based on what type of model for true_beta.
  # if no type given, assume multiple linear regression
  if(!is.null(seed)) set.seed(seed)
  if(is.null(type)){ # multidimensional case
    if(dim(D)[2] != length(true_beta)) stop("dimension of D doesn't match length of true_beta in case where type == NULL")
    X = D
  } else{ # regress on transformations of x
    X = constructDesignX(D, N, type)
  }
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    # Y[ , j] = rmvnorm(1, X %*% true_beta, var_e * diag(N))
    Y[ , j] = X %*% true_beta + sqrt(var_e) * rnorm(N)
  }
  return(Y)
}

simulateYvs = function(D, N, true_beta, var_e, numSims, seed = NULL){
  simulateY(D, N, true_beta, var_e, numSims, type = NULL, seed = NULL)
}

simulateYN = function(X, true_beta, var_e, numSims, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  if(dim(X)[2] != length(true_beta)) stop("simulateYN : X and true_beta aren't compatible!")
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  # N x numSims -- each row is a simulation of N data points, 
  # where N is the number of rows in design matrix X
  Y = rmvnorm(numSims, X %*% true_beta, var_e * diag(N))
}

simulateY_fromfunction = function(x, f, var_e, numSims, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  Y = f(x) + rnorm(length(x), 0, sqrt(var_e))
}