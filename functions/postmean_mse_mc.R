# require("construct_design_matrix.R")
# require("simulate_y.R")

# first compute estimator, posterior mean
getPostMean = function(y, D, N, true_beta, var_e, var_mean, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  D_postvar = postvar(D, N, var_e, var_mean, type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(var_mean, true_beta))
  return(D_postmean)
}

getEmpMSE = function(postmean_beta, true_beta){
  return((postmean_beta - true_beta)^2)
}

calcExpPostMeanMSE_old = function(D, N, true_beta, numSims, mean_beta, var_e, var_mean, type, diagPrior = TRUE){
  Ysims = simulateY(D, N, mean_beta, var_mean, var_e, numSims, type = type)
  post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, true_beta, var_e, var_mean, type, diagPrior = diagPrior))
  empMSEs = apply(post_means, 2, FUN = function(x) getEmpMSE(x, true_beta))
  expEpiricalMSE = apply(empMSEs, 1, mean)
  return(expEpiricalMSE)
}

calcExpPostMeanMSE = function(D, N, true_beta, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                              numSims = 100, true_type = NULL, type = NULL, seed = 123){
  set.seed(seed)
  Ysims = simulateY(D, N, true_beta, var_e, numSims, true_type, seed)
  post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, true_beta, var_e, var_mean, type, diagPrior = diagPrior))
  empMSEs = apply(post_means, 2, FUN = function(x) getEmpMSE(x, true_beta))
  expEpiricalMSE = apply(empMSEs, 1, mean)
  return(expEpiricalMSE)
}