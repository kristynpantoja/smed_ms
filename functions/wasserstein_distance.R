# Wasserstein distance betwen two normals
Wasserstein_distance = function(mu1, mu2, var1, var2, dim = 1){
  # Normal(mu1, var1)
  # Normal(mu2, var2)
  # dim:
  #   1 is for univariate case
  #   > 1 is for multivariate case
  if(dim == 1){
    wass = sqrt((mu1 - mu2)^2 + (sqrt(var1) - sqrt(var2))^2)
  } else{
    if(dim > 1){
      sqrt_var2 = sqrtm(var2)
      wass = sqrt(crossprod(mu1 - mu2) + sum(diag(var1 + var2 - 2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2))))
    } else{
      stop("invalid dim")
    }
  }
  return(as.numeric(wass))
}



# require("posterior_mean.R")
Wasserstein_distance_postpred = function(x, initD, y, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                         f0, f1, type, N2, dim = 1){
  initN = length(initD)
  
  # posterior distribution of beta
  x0 = constructDesignX(x, 1, type[1])
  X0 = constructDesignX(initD, initN, type[1])
  postvar0 = postvar(initD, initN, var_e, var_beta0, type[1])
  postmean0 = postmean(y, initD, initN, mean_beta0, var_beta0, var_e, type[1])
  
  x1 = constructDesignX(x, 1, type[2])
  X1 = constructDesignX(initD, initN, type[2])
  postvar1 = postvar(initD, initN, var_e, var_beta1, type[2])
  postmean1 = postmean(y, initD, initN, mean_beta1, var_beta1, var_e, type[2])
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(x0) %*% postmean0
  postpredvar0 = x0 %*% postvar0 %*% t(x0) + var_e
  
  postpredmu1 = t(x1) %*% postmean1
  postpredvar1 = x1 %*% postvar1 %*% t(x1) + var_e
  
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1)
  return(as.numeric(wass))
}