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
      stop("error in WN() : invalid dim")
    }
  }
  return(as.numeric(W))
}


####################################################################
### SeqMED (uses posterior predictive distribution of y) ###########
####################################################################

# 1 dimension, or with transformations of x
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

# multidimensional, for variable selection
WNlmvs = function(
  x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var
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

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# compare covariance functions or other parameters in 1 dimension
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

WNgpvs = function(x, Kinv0, Kinv1, initD0, initD1, y, model0, model1){
  x = t(as.matrix(x))
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