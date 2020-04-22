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

# X = constructDesignX(x_random2[ , indicesT], Ntot2, NULL)
# y1 = as.vector(rmvnorm(1, X %*% betaT, sigmasq * diag(Ntot2)))
# y2 = X %*% betaT + sqrt(sigmasq) * rnorm(Ntot2, 0, 1)
# y3 = X %*% betaT + rnorm(Ntot2, 0, sqrt(sigmasq))
# lm(y1 ~ -1 + X)
# lm(y2 ~ -1 + X)
# lm(y3 ~ -1 + X)
# microbenchmark(
#   as.vector(rmvnorm(1, X %*% betaT, sigmasq * diag(Ntot2))),
#   X %*% betaT + sqrt(sigmasq) * rnorm(Ntot2, 0, 1),
#   X %*% betaT + rnorm(Ntot2, 0, sqrt(sigmasq))
# )
