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


####################################################################
### SeqMED (uses posterior predictive distribution of y) ###########
####################################################################
# require("posterior_mean.R")
# require("construct_design_matrix.R")
# require("posterior_variance.R")

# 1 dimension, or with transformations of x
# formerly named Wasserstein_distance_postpred
WNlm = function(
  x, postmean0, postmean1, postvar0, postvar1, model0, model1, error.var, 
  dim = 1
){
  x0 = t(model0$designMat(x))
  x1 = t(model1$designMat(x))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(x0) %*% postmean0
  postpredy_var0 = t(x0) %*% postvar0 %*% x0 + error.var
  
  postpredy_mu1 = t(x1) %*% postmean1
  postpredy_var1 = t(x1) %*% postvar1 %*% x1 + error.var
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}
WNlm_old = function(
  x, postmean0, postmean1, postvar0, postvar1, error.var, type, dim = 1){
  x0 = t(constructDesignX(x, 1, type[1]))
  x1 = t(constructDesignX(x, 1, type[2]))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(x0) %*% postmean0
  postpredy_var0 = t(x0) %*% postvar0 %*% x0 + error.var
  
  postpredy_mu1 = t(x1) %*% postmean1
  postpredy_var1 = t(x1) %*% postvar1 %*% x1 + error.var
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

# multidimensional, for variable selection
# formerly named Wasserstein_vs
WNlmvs = function(
  x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var, 
  dim = 1
){
  x = t(as.matrix(x))
  
  x0 = x[ , model0$indices, drop = FALSE]
  x1 = x[ , model1$indices, drop = FALSE]
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = x0 %*% postmean0
  postpredy_var0 = x0 %*% postvar0 %*% t(x0) + error.var
  
  postpredy_mu1 = x1 %*% postmean1
  postpredy_var1 = x1 %*% postvar1 %*% t(x1) + error.var
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}
WNlmvs_old = function(
  x, variables0, variables1, postmean0, postmean1, postvar0, postvar1, signal.var, 
  dim = 1){
  x0 = t(x)[ , variables0]
  x1 = t(x)[ , variables1]
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(x0) %*% postmean0
  postpredy_var0 = t(x0) %*% postvar0 %*% x0 + signal.var
  
  postpredy_mu1 = t(x1) %*% postmean1
  postpredy_var1 = t(x1) %*% postvar1 %*% x1 + signal.var
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# compare covariance functions or other parameters in 1 dimension
# formerly named Wasserstein_distance_postpred_gp
WNgp = function(x, Kinv0, Kinv1, initD, y, model0, model1){
  
  # posterior distribution of beta
  k0 = t(as.matrix(getCov(
    X1 = x, X2 = initD, type = model0$type, l = model0$l, p = model0$p, 
    signal.var = model0$signal.var)))
  k1 = t(as.matrix(getCov(
    X1 = x, X2 = initD, type = model1$type, l = model1$l, p = model1$p, 
    signal.var = model1$signal.var)))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(k0) %*% Kinv0 %*% y
  postpredy_var0 = 1 - t(k0) %*% Kinv0 %*% k0
  if(postpredy_var0 < 0) postpredy_var0 = 0 # only happens when too-small
  
  postpredy_mu1 = t(k1) %*% Kinv1 %*% y
  postpredy_var1 = 1 - t(k1) %*% Kinv1 %*% k1
  if(postpredy_var1 < 0) postpredy_var1 = 0 # same reason
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}

# really slow, because have to re-compute Kinv0, Kinv1
WNgp.new = function(x, Kinv0, Kinv1, initD, y, model0, model1){
  
  gp0 = getGPPredictive(
    x, initD, y, model0$type, model0$l, model0$p, model0$signal.var, 
    model0$error.var)
  gp1 = getGPPredictive(
    x, initD, y, model1$type, model1$l, model1$p, model1$signal.var, 
    model1$error.var)
  
  # error checking
  if(gp0$pred_var < 0) gp0$pred_var = 0 # only happens when too-small
  if(gp1$pred_var < 0) gp1$pred_var = 0
  
  W = WN(gp0$pred_mean, gp0$pred_mean, gp0$pred_var, gp1$pred_var)
  return(as.numeric(W))
}

# multidimensional, for variable selection
# formerly named Wasserstein_distance_postpred_gpvs
# WNgpvs = function(
#   x, Kinv0, Kinv1, variables0, variables1, initD0, initD1, y, signal.var, type, l){
#   x = t(as.matrix(x))
#   
#   # posterior distribution of beta
#   k0 = t(as.matrix(getCov(x[ , variables0, drop = FALSE], 
#                           initD0, type[1], l[1])))
#   k1 = t(as.matrix(getCov(x[ , variables1, drop = FALSE], 
#                           initD1, type[2], l[2])))
#   
#   # posterior predictive distribution of y, for candidate x
#   postpredy_mu0 = t(k0) %*% Kinv0 %*% y
#   postpredy_var0 = signal.var * (1 - t(k0) %*% Kinv0 %*% k0)
#   if(postpredy_var0 < 0) postpredy_var0 = 0 
#   
#   postpredy_mu1 = t(k1) %*% Kinv1 %*% y
#   postpredy_var1 = signal.var * (1 - t(k1) %*% Kinv1 %*% k1)
#   if(postpredy_var1 < 0) postpredy_var1 = 0 
#   W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1, dim = 1)
#   return(as.numeric(W))
# }

WNgpvs = function(x, Kinv0, Kinv1, initD0, initD1, y, model0, model1){
  x = t(as.matrix(x))
  # print(x)
  # posterior distribution of beta
  k0 = t(as.matrix(getCov(
    X1 = x[, model0$indices, drop = FALSE], X2 = initD0, type = model0$type, 
    l = model0$l, p = model0$p, signal.var = model0$signal.var)))
  k1 = t(as.matrix(getCov(
    X1 = x[, model1$indices, drop = FALSE], X2 = initD1, type = model1$type, 
    l = model1$l, p = model1$p, signal.var = model1$signal.var)))
  
  # posterior predictive distribution of y, for candidate x
  postpredy_mu0 = t(k0) %*% Kinv0 %*% y
  postpredy_var0 = 1 - t(k0) %*% Kinv0 %*% k0
  if(postpredy_var0 < 0) postpredy_var0 = 0 # only happens when too-small
  
  postpredy_mu1 = t(k1) %*% Kinv1 %*% y
  postpredy_var1 = 1 - t(k1) %*% Kinv1 %*% k1
  if(postpredy_var1 < 0) postpredy_var1 = 0 # same reason
  
  W = WN(postpredy_mu0, postpredy_mu1, postpredy_var0, postpredy_var1)
  return(as.numeric(W))
}
