# require("construct_design_matrix.R")

simulateY = function(D, N, true_beta, var_e, numSims, type = NULL, seed = NULL){
  # note that type is based on what type of model for true_beta.
  # if no type given, assume multiple linear regression
  if(!is.null(seed)) set.seed(seed)
  if(is.null(type)){
    if(dim(D)[2] != length(true_beta)) stop("dimension of D doesn't match length of true_beta in case where type == NULL")
    X = D
  } else{
    X = constructDesignX(D, N, type)
  }
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    Y[ , j] = rmvnorm(1, X %*% true_beta, var_e * diag(N))
  }
  return(Y)
}
