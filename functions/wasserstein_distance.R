# Wasserstein distance betwen two normals
# formerly named Wasserstein_distance
WN = function(mu1, mu2, var1, var2, dim = 1){
  # Normal(mu1, var1)
  # Normal(mu2, var2)
  # dim:
  #   1 is for univariate case
  #   > 1 is for multivariate case
  if(dim == 1){
    W = sqrt((mu1 - mu2)^2 + (sqrt(var1) - sqrt(var2))^2)
  } else{
    if(dim > 1){
      sqrt_var2 = sqrtm(var2)
      W = sqrt(crossprod(mu1 - mu2) + 
                 sum(diag(var1 + var2 - 
                            2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2)
                          )
                     )
               )
    } else{
      stop("error in WassN : invalid dim")
    }
  }
  return(as.numeric(W))
}



# require("posterior_mean.R")
# require("construct_design_matrix.R")
# require("posterior_variance.R")

####################################################################
### SeqMED (uses posterior predictive distribution of y) ###########
####################################################################

# 1 dimension, or with transformations of x
# formerly named Wasserstein_distance_postpred
WNlm = function(
  x, postmean0, postmean1, postvar0, postvar1, var_e, type, dim = 1){
  x0 = t(constructDesignX(x, 1, type[1]))
  x1 = t(constructDesignX(x, 1, type[2]))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(x0) %*% postmean0
  postpredy_var0 = t(x0) %*% postvar0 %*% x0 + var_e
  
  postpredy_mu1 = t(x1) %*% postmean1
  postpredy_var1 = t(x1) %*% postvar1 %*% x1 + var_e
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

# multidimensional, for variable selection
# formerly named Wasserstein_vs
WNlmvs = function(
  x, variables0, variables1, postmean0, postmean1, postvar0, postvar1, var_e, 
  dim = 1){
  x0 = t(constructDesignX(x, 1, NULL))[ , variables0]
  x1 = t(constructDesignX(x, 1, NULL))[ , variables1]
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(x0) %*% postmean0
  postpredy_var0 = t(x0) %*% postvar0 %*% x0 + var_e
  
  postpredy_mu1 = t(x1) %*% postmean1
  postpredy_var1 = t(x1) %*% postvar1 %*% x1 + var_e
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# compare covariance functions or other parameters in 1 dimension
# formerly named Wasserstein_distance_postpred_gp
WNgp = function(x, Kinv0, Kinv1, initD, y, var_e, type, l){
  
  # posterior distribution of beta
  k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
  k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(k0) %*% Kinv0 %*% y
  postpredy_var0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
  if(postpredy_var0 < 0) postpredy_var0 = 0 # !!!!!!!!!! #############################################
  
  postpredy_mu1 = t(k1) %*% Kinv1 %*% y
  postpredy_var1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
  if(postpredy_var1 < 0) postpredy_var1 = 0 # !!!!!!!!!! #############################################
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

# multidimensional, for variable selection
# formerly named Wasserstein_distance_postpred_gpvs
WNgpvs = function(
  x, Kinv0, Kinv1, variables0, variables1, initD0, initD1, y, var_e, type, l){
  x = t(as.matrix(x))
  
  # posterior distribution of beta
  k0 = t(as.matrix(getCov(x[ , variables0, drop = FALSE], 
                          initD0, type[1], l[1])))
  k1 = t(as.matrix(getCov(x[ , variables1, drop = FALSE], 
                          initD1, type[2], l[2])))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(k0) %*% Kinv0 %*% y
  postpredy_var0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
  if(postpredy_var0 < 0) postpredy_var0 = 0 # !!!!!!!!!! ##########################################
  
  postpredy_mu1 = t(k1) %*% Kinv1 %*% y
  postpredy_var1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
  if(postpredy_var1 < 0) postpredy_var1 = 0 # !!!!!!!!!! ############################################
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1, dim = 1)
  return(as.numeric(W))
}