require("wasserstein_distance.R")
require("charge_function_q.R")

# 1d Covariance functions
# radial aka squared exponential aka gaussian
C_fn_elementwise_sqexp1d = function(Xi,Xj, l) exp(-0.5 * ((Xi - Xj) / l)^2)
# exponential, a.k.a. matern, v = 1/2
C_fn_elementwise_exp1d = function(Xi, Xj, l) exp(- abs(Xi - Xj) / l)
# matern, v = 3/2
C_fn_elementwise_matern1dv1 = function(Xi, Xj, l){
  r = abs(Xi - Xj)
  (1 + sqrt(5) * r / l + 5 * r^2 / (3 * l^2)) * exp(- sqrt(3) * r / l)
}
# compute covariance matrix
C_fn_1d = function(X, Y, l, C_fn_type = 1){
  C_fn_elementwise = NULL
  if(C_fn_type == 1) C_fn_elementwise = C_fn_elementwise_sqexp1d
  else if(C_fn_type == 2) C_fn_elementwise = C_fn_elementwise_exp1d
  else if(C_fn_type == 3 ) C_fn_elementwise = C_fn_elementwise_matern1dv1
  else stop("invalid k argument, i.e. type specification of covariance function for GP")
  outer(X, Y, FUN = C_fn_elementwise, l)
}

getPredDistr = function(x_star, x_train, y_train, l, C_fn_type){
  k_star = t(as.matrix(C_fn_1d(x_star, x_train, l, C_fn_type)))
  K_obs = C_fn_1d(x_train, x_train, l, C_fn_type)
  pred_mean = t(k_star) %*% solve(K_obs, y_train)
  pred_cov = C_fn_1d(x_star, x_star, l, C_fn_type) - t(k_star) %*% solve(K_obs, k_star)
  return(list("pred_mean" = pred_mean, "pred_cov" = pred_cov))
}

getW = function(x_star, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1){
  pred_distr0 = getPredDistr(x_star, x_train, y_train, l0, C_fn_type0)
  pred_distr1 = getPredDistr(x_star, x_train, y_train, l1, C_fn_type1)
  mu1 = pred_distr0$pred_mean # mean of marginal dist of y | H1
  mu2 = pred_distr1$pred_mean
  var1 = pred_distr1$pred_cov
  var2 = pred_distr1$pred_cov
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2, dim = 1)
  return(Wass_dist)
}



###################
### one at time ###
###################

f_min_gp = function(candidate, D, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                 f0, f1, type, var_margy0, var_margy1, p, log_space = FALSE){
  if(log_space == FALSE) {
    result = q_gp(candidate, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1)^k * sum(sapply(D, function(x_i) (q(x_i, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1) / sqrt((x_i - candidate)^2))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q(candidate, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1)) + k * log(q(x_i, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1)) - k * log(sqrt((x_i - candidate)^2)))
    result = exp(logSumExp(terms))
    return(result)
  }
}

# one at a time algorithm
MED_ms_gp_1pt = function(x_train, y_train, l0, l1, C_fn_type0, C_fn_type1, 
                   N = 11, numCandidates = 10^5, k = 4, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
                   genCandidates = 1, initialpt = 2){
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  # remove candidate points that are already in the design
  # this is necessary bc they lead to a posterior predictive distribution with
  # standard deviation/covariance = 0 and means equal to a design point
  # -> W = 0 -> division by 0 in q evaluation
  candidates = candidates[!(candidates %in% x_train)]
  # but also can't have candidate points that are too close to them either...
  # lead to negative variances in predictive distribution somehow
  
  # -- Initialize 1st Design Point in D -- #
  w_evals = sapply(candidates, FUN = function(x) getW(x, x_train, y_train, l0, l1, C_fn_type0, C_fn_type1))
  xinitind = which.max(w_evals)
  
  
}

