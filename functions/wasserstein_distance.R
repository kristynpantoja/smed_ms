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
# require("construct_design_matrix.R")
# require("posterior_variance.R")

####################################################################
### MMED with data (uses posterior predictive distribution of y) ###
####################################################################

# 1 dimension, or with transformations of x
Wasserstein_distance_postpred = function(x, postmean0, postmean1, postvar0, postvar1, var_e, type, dim = 1){
  
  # posterior distribution of beta
  x0 = t(constructDesignX(x, 1, type[1]))
  x1 = t(constructDesignX(x, 1, type[2]))
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(x0) %*% postmean0
  postpredvar0 = t(x0) %*% postvar0 %*% x0 + var_e
  
  postpredmu1 = t(x1) %*% postmean1
  postpredvar1 = t(x1) %*% postvar1 %*% x1 + var_e
  
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1)
  return(as.numeric(wass))
}

# multidimensional 
Wasserstein_distance_postpred_multidim = function(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, dim = 1){
  # Wasserstein_distance_postpred(x, postmean0, postmean1, postvar0, postvar1, var_e, type = c(NULL, NULL), dim = 1)
  
  # posterior distribution of beta
  x0 = t(constructDesignX(x, 1, NULL))[ , indices0]
  x1 = t(constructDesignX(x, 1, NULL))[ , indices1]
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(x0) %*% postmean0
  postpredvar0 = t(x0) %*% postvar0 %*% x0 + var_e
  
  postpredmu1 = t(x1) %*% postmean1
  postpredvar1 = t(x1) %*% postvar1 %*% x1 + var_e
  
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1)
  return(as.numeric(wass))
}

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# compare covariance functions or other parameters in 1 dimension
Wasserstein_distance_postpred_gp = function(x, Kinv0, Kinv1, initD, y, var_e, type, l){
  
  # posterior distribution of beta
  k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
  k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(k0) %*% Kinv0 %*% y
  postpredvar0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
  if(postpredvar0 < 0) postpredvar0 = 0 # !!!!!!!!!!
  
  postpredmu1 = t(k1) %*% Kinv1 %*% y
  postpredvar1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
  if(postpredvar1 < 0) postpredvar1 = 0 # !!!!!!!!!!
  
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1)
  return(as.numeric(wass))
}

################################################
### MMED GP for Variable Selection, 1D vs 2D ###
################################################

# meant to be able to handle 2d dimensional input, for variable selection problem
# Wasserstein_distance_postpred_gpvsold = function(x, Kinv0, Kinv1, subdim, subinitD = NULL, initD, y, var_e, type, l){
#   x = as.matrix(x)
#   if(dim(x)[1] !=1) x = t(x)
#   
#   # posterior distribution of beta
#   if(is.null(subinitD)){
#     k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
#   } else{
#     # k0 = t(as.matrix(getCov(x[ , subdim, drop = FALSE], subinitD, type[1], l[1])))
#     k0 = t(as.matrix(getCov(x, subinitD, type[1], l[1])))
#   }
#   k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
#   
#   # posterior predictive distribution of y, for candidate x
#   postpredmu0 = t(k0) %*% Kinv0 %*% y
#   postpredvar0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
#   if(postpredvar0 < 0) postpredvar0 = 0 # !!!!!!!!!!
#   
#   postpredmu1 = t(k1) %*% Kinv1 %*% y
#   postpredvar1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
#   if(postpredvar1 < 0) postpredvar1 = 0 # !!!!!!!!!!
#   wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1, dim = 1)
#   return(as.numeric(wass))
# }
Wasserstein_distance_postpred_gpvs = function(x, Kinv0, Kinv1, subdim, subinitD = NULL, initD, y, var_e, type, l){
  if(!is.null(subinitD)) if(typeof(subinitD) == "list") subinitD = as.matrix(subinitD)
  if(typeof(initD) == "list") initD = as.matrix(initD)
  x = as.matrix(x)
  if(dim(x)[1] !=1) x = t(x)
  
  # posterior distribution of beta
  if(is.null(subinitD)){
    k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
  } else{
    k0 = t(as.matrix(getCov(x[ , subdim, drop = FALSE], subinitD, type[1], l[1])))
  }
  k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(k0) %*% Kinv0 %*% y
  postpredvar0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
  if(postpredvar0 < 0) postpredvar0 = 0 # !!!!!!!!!!
  
  postpredmu1 = t(k1) %*% Kinv1 %*% y
  postpredvar1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
  if(postpredvar1 < 0) postpredvar1 = 0 # !!!!!!!!!!
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1, dim = 1)
  return(as.numeric(wass))
}