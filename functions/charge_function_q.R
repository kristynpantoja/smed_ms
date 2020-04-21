# charge function evaluated at x

##########
### 1D ###
##########

q_nodata = function(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, var_margy0, var_margy1, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 1
  mu1.temp = f0(x) # mean of marginal dist of y | H0
  mu2.temp = f1(x) # mean of marginal dist of y | H1
  var1.temp = var_marginaly(x, var_beta0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2.temp = var_marginaly(x, var_beta1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1.temp, mu2.temp, var1.temp, var2.temp)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

q_data = function(x, postmean0, postmean1, postvar0, postvar1, var_e, type, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 1
  
  Wass_dist = Wasserstein_distance_postpred(x, postmean0, postmean1, postvar0, postvar1, var_e, type)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

q_vs = function(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, alpha = NULL, buffer = 0){
  if(is.null(alpha)) alpha = 1
  
  Wass_dist = Wasserstein_vs(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

# q_data_multidim = function(x, postmean0, postmean1, postvar0, postvar1, var_e, type, p, alpha = NULL, buffer = 0){
#   if(length(type) != 2) stop("type should be vector with length == 2")
#   if(is.null(alpha)) alpha = 1
#   Wass_dist = Wasserstein_distance_postpred_multidim(x, postmean0, postmean1, postvar0, postvar1, var_e, type)
#   q_exponent = alpha / (2 * p)
#   return(1.0 / (Wass_dist + buffer)^q_exponent)
# }

q_gp = function(x, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 1
  Wass_dist = Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

# q_gpvs = function(x, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l, p, alpha = NULL, buffer = 0){
#   if(length(type) != 2) stop("type should be vector with length == 2")
#   if(is.null(alpha)) alpha = 1
#   Wass_dist = Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l)
#   q_exponent = alpha / (2 * p)
#   return(1.0 / (Wass_dist + buffer)^q_exponent)
# }
q_gpvs = function(x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, var_e, type, l, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 1
  Wass_dist = Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, var_e, type, l)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

##########
### 2D ###
##########

# charge function evaluated at x
q_nodata_2d = function(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, var_margy0, var_margy1, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 1
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly_2d(x, var_beta0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly_2d(x, var_beta1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}