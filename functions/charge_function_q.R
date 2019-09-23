# charge function evaluated at x

##########
### 1D ###
##########

q = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 2 * p
  mu1.temp = f0(x) # mean of marginal dist of y | H0
  mu2.temp = f1(x) # mean of marginal dist of y | H1
  var1.temp = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2.temp = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1.temp, mu2.temp, var1.temp, var2.temp)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

##########
### 2D ###
##########

# charge function evaluated at x
q_2d = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, alpha = NULL, buffer = 0){
  if(length(type) != 2) stop("type should be vector with length == 2")
  if(is.null(alpha)) alpha = 2 * p
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly_2d(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly_2d(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  q_exponent = alpha / (2 * p)
  return(1.0 / (Wass_dist + buffer)^q_exponent)
}

##########
### GP ###
##########

# charge function evaluated at x
q_gp = function(x_star, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1){
  pred_distr0 = getPredDistr(x_star, x_train, y_train, l0, C_fn_type0)
  pred_distr1 = getPredDistr(x_star, x_train, y_train, l1, C_fn_type1)
  mu1 = pred_distr0$pred_mean # mean of marginal dist of y | H1
  mu2 = pred_distr1$pred_mean
  var1 = pred_distr0$pred_cov
  var2 = pred_distr1$pred_cov
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2, dim = 1)
  return(1.0 / Wass_dist)
}
