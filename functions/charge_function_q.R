# charge function evaluated at x

# for D0 in normal polynomial models with one explanatory variable, 
#   when there is no data
q_mmed = function(
  x, mean_beta0, mean_beta1, var_beta0, var_beta1, error.var, f0, f1, type, 
  var_margy0, var_margy1, p, alpha = 1
){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1.temp = f0(x) # mean of marginal dist of y | H0
  mu2.temp = f1(x) # mean of marginal dist of y | H1
  var1.temp = var_marginaly(x, var_beta0, error.var, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2.temp = var_marginaly(x, var_beta1, error.var, type = type[2], var_margy1) # variance of marginal dist of y | H1
  W = WN(mu1.temp, mu2.temp, var1.temp, var2.temp)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}

# for D1, ..., DT in normal polynomial models with one explanatory variable
q_seqmed = function(
  x, postmean0, postmean1, postvar0, postvar1, error.var, type, p = 1, alpha = 1
){
  if(length(type) != 2) stop("type should be vector with length == 2")
  W = WNlm(x, postmean0, postmean1, postvar0, postvar1, error.var, type)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}

# for D1, ..., DT in normal multiple linear regression models
q_vs = function(
  x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, error.var, p, 
  alpha = 1
){
  W = WNlmvs(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, error.var)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}

# for D1, ..., DT in gp models with one explanatory variable
q_gp = function(
  x, Kinv0, Kinv1, initD, y, signal.var, type, l, p, alpha = 1, buffer = 0
){
  if(length(type) != 2) stop("type should be vector with length == 2")
  W = WNgp(x, Kinv0, Kinv1, initD, y, signal.var, type, l) + buffer
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}

# for D1, ..., DT in multivariate gp models
q_gpvs = function(
  x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, type, l, p, 
  alpha = 1
){
  if(length(type) != 2) stop("type should be vector with length == 2")
  W = WNgpvs(x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, type, l)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}
