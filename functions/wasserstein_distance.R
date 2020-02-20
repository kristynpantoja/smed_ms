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


# Wasserstein_distance_gp = function(mu1, mu2, var1, var2, dim = 1){
#   # Normal(mu1, var1)
#   # Normal(mu2, var2)
#   # dim:
#   #   1 is for univariate case
#   #   > 1 is for multivariate case
#   if(dim == 1){
#     wass = tryCatch({sqrt((mu1 - mu2)^2 + (sqrt(var1) - sqrt(var2))^2)}, 
#                     warning = function(warn) {0})
#   } else{
#     if(dim > 1){
#       sqrt_var2 = sqrtm(var2)
#       wass = tryCatch({sqrt(crossprod(mu1 - mu2) + 
#                               sum(diag(var1 + var2 - 
#                                          2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2))))}, 
#                       error = function(e) {0})
#         
#     } else{
#       stop("invalid dim")
#     }
#   }
#   return(as.numeric(wass))
# }


# require("posterior_mean.R")
# require("construct_design_matrix.R")
# require("posterior_variance.R")


Wasserstein_distance_postpred2 = function(x, postmean0, postmean1, postvar0, postvar1, var_e, type, dim = 1){
  
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

Wasserstein_distance_postpred_gp2 = function(x, Kinv0, Kinv1, subinitD = NULL, initD, y, var_e, type, l){
  x = as.matrix(x)
  if(dim(x)[1] !=1) x = t(x)
  # print(x)
  
  # posterior distribution of beta
  if(is.null(subinitD)){
    k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
  } else{
    k0 = t(as.matrix(getCov(x, subinitD, type[1], l[1])))
  }
  k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
  
  # posterior predictive distribution of y, for candidate x
  postpredmu0 = t(k0) %*% Kinv0 %*% y
  postpredvar0 = var_e * (1 - t(k0) %*% Kinv0 %*% k0)
  if(postpredvar0 < 0) postpredvar0 = 0 # !!!!!!!!!!
  
  postpredmu1 = t(k1) %*% Kinv1 %*% y
  postpredvar1 = var_e * (1 - t(k1) %*% Kinv1 %*% k1)
  if(postpredvar1 < 0) postpredvar1 = 0 # !!!!!!!!!!
  # print(paste("dim of mu1: ", paste(dim(postpredmu0))))
  # print(paste("dim of mu2: ", dim(postpredmu1)))
  # print(paste("dim of var1: ", dim(postpredvar0)))
  # print(paste("dim of var1: ", dim(postpredvar1)))
  # print(paste("mu1: ", postpredmu0))
  # print(paste("mu2: ", postpredmu1))
  # print(paste("var1: ", postpredvar0))
  # print(paste("var1: ", postpredvar1))
  wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1, dim = 1)
  return(as.numeric(wass))
}


# Wasserstein_distance_postpred_gp_NEW = function(x, Kc0, Kc1, initD, y, var_e, type, l, dim = 1){
#   
#   # posterior distribution of beta
#   k0 = t(as.matrix(getCov(x, initD, type[1], l[1])))
#   k1 = t(as.matrix(getCov(x, initD, type[2], l[2])))
#   
#   # posterior predictive distribution of y, for candidate x
#   postpredmu0 = t(k0) %*% forwardsolve(Kc0, y)
#   postpredvar0 = var_e * (1 - t(k0) %*% forwardsolve(Kc0, k0))
#   if(postpredvar0 < 0) postpredvar0 = 0 # !!!!!!!!!!
#   
#   postpredmu1 = t(k1) %*% forwardsolve(Kc1, y)
#   postpredvar1 = var_e * (1 - t(k1) %*% forwardsolve(Kc1, k1))
#   if(postpredvar1 < 0) postpredvar1 = 0 # !!!!!!!!!!
#   
#   wass = Wasserstein_distance(postpredmu0, postpredmu1, postpredvar0, postpredvar1)
#   return(as.numeric(wass))
# }
