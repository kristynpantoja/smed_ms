# charge function evaluated at x

# for D0 in normal polynomial models with one explanatory variable, 
#   when there is no data
q_mmed = function(
  x, model0, model1, error.var, var_margy0, var_margy1, p, alpha = 1
){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1.temp = f0(x) # mean of marginal dist of y | H0
  mu2.temp = f1(x) # mean of marginal dist of y | H1
  var1.temp = var_margy0(x) # variance of marginal dist of y | H0
  var2.temp = var_margy1(y) # variance of marginal dist of y | H1
  W = WN(mu1.temp, mu2.temp, var1.temp, var2.temp)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
}

# for D1, ..., DT in normal polynomial models with one explanatory variable
q_seqmed = function(
  x, postmean0, postmean1, postvar0, postvar1, model0, model1, error.var, 
  p = 1, alpha = 1
){
  
  W = WNlm(
    x = x, postmean0 = postmean0, postmean1 = postmean1, 
    postvar0 = postvar0, postvar1 = postvar1, 
    model0 = model0, model1 = model1, error.var = error.var)
  q_exponent = alpha / (2 * p)
  return(1.0 / W^q_exponent)
}

# for D1, ..., DT in normal multiple linear regression models
q_vs = function(
  x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var, p, 
  alpha = 1
){
  W = WNlmvs(x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var)
  q_exponent = alpha / (2 * p)
  return(1.0 / (W)^q_exponent)
} 

# for D1, ..., DT in gp models with one explanatory variable
q_gp = function(
  x, Kinv0, Kinv1, initD, y, p, alpha = 1, buffer = 0, model0, model1
){
  W = WNgp(x, Kinv0, Kinv1, initD, y, model0, model1) + buffer
  q_exponent = alpha / (2 * p)
  return(1.0 / W^q_exponent)
}
qcap_gp = function(
  x, Kinv0, Kinv1, initD, y, p, alpha = 1, model0, model1, cap = 1e20
){
  q_val = q_gp(
    x = x, Kinv0 = Kinv0, Kinv1 = Kinv1, initD = initD, y = y, p = p, 
    alpha = alpha, model0 = model0, model1 = model1)
  q_prime = min(q_val, cap)
  return(q_prime)
}

# # for D1, ..., DT in multivariate gp models
q_gpvs = function(
  x, Kinv0, Kinv1, initD0, initD1, y, p, alpha = 1, buffer = 0, model0, model1
){
  W = WNgpvs(x, Kinv0, Kinv1, initD0, initD1, y, model0, model1) + buffer
  q_exponent = alpha / (2 * p)
  return(1.0 / W^q_exponent)
}
qcap_gpvs = function(
  x, Kinv0, Kinv1, initD0, initD1, y, p, alpha = 1, model0, model1, cap = 1e20
){
  q_val = q_gpvs(
    x = x, Kinv0 = Kinv0, Kinv1 = Kinv1, initD0 = initD0, initD1 = initD1, 
    y = y, p = p, alpha = alpha, model0 = model0, model1 = model1)
  q_prime = min(q_val, cap)
  return(q_prime)
}
