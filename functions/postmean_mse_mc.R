# first compute estimator, posterior mean
getPostMean = function(y, D, N, mean_beta, var_e, var_mean, type, diagPrior = TRUE){
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  D_postvar = postvar(D, N, var_e, var_mean, type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(var_mean, mean_beta))
  return(D_postmean)
}

getEmpMSE = function(postmean, truemean){
  return((postmean - truemean)^2)
}

calcExpPostMeanMSE = function(D, N, true_beta, numSims, mean_beta, var_e, var_mean, type, diagPrior = TRUE){
  Ysims = simulateY(D, N, mean_beta, var_mean, var_e, numSims, type = type)
  post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D = D, N = N, mean_beta = mean_beta, var_e = var_e, var_mean = var_mean, type = type, diagPrior = diagPrior))
  empMSEs = apply(post_means, 2, FUN = function(x) getEmpMSE(x, true_beta))
  expEpiricalMSE = apply(empMSEs, 1, mean)
  return(expEpiricalMSE)
}
