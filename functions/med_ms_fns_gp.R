# require("wasserstein_distance.R")
# require("charge_function_q.R")

# 1d Covariance functions
# radial aka squared exponential aka gaussian
phi_sqexp = function(Xi,Xj, l) exp(-0.5 * ((Xi - Xj) / l)^2)
# exponential, a.k.a. matern, v = 1/2
phi_exp = function(Xi, Xj, l) exp(- abs(Xi - Xj) / l)
# matern, v = 3/2
phi_matern2 = function(Xi, Xj, l){
  r = abs(Xi - Xj)
  (1 + sqrt(5) * r / l + 5 * r^2 / (3 * l^2)) * exp(- sqrt(3) * r / l)
}
# general matern, l = 1
phi_matern = function(Xi, Xj, l, nu = 1) RandomFieldsUtils::matern(abs(Xi - Xj), nu, scaling = "matern")

# compute covariance matrix
getCov = function(X, Y, type = 1, l = 1){ # change l = 1 to a list theta = NULL
  phi = NULL
  if(type == 1) phi = phi_sqexp
  else if(type == 2) phi = phi_matern
  else if(type == 3) phi =phi_exp
  else if(type == 4) phi = phi_matern2
  else stop("invalid type specification of covariance function for GP")
  outer(X, Y, FUN = phi, l)
}

getPredDistr = function(x_star, x_train, y_train, type, l, nugget){
  k_star = t(as.matrix(getCov(x_star, x_train, type, l)))
  if(is.null(nugget)) K_obs = getCov(x_train, x_train, type, l)
  else K_obs = getCov(x_train, x_train, type, l) + diag(rep(nugget, length(x_train)))
  pred_mean = t(k_star) %*% solve(K_obs, y_train)
  pred_cov = getCov(x_star, x_star, type, l) - t(k_star) %*% solve(K_obs, k_star)
  return(list("pred_mean" = as.vector(pred_mean), "pred_var" = as.vector(pred_cov)))
}

getW = function(x_star, x_train, y_train, type01, l01, nugget){
  pred_distr0 = getPredDistr(x_star, x_train, y_train, type01[1], l01[1], nugget)
  pred_distr1 = getPredDistr(x_star, x_train, y_train, type01[2], l01[2], nugget)
  mu1 = pred_distr0$pred_mean # mean of marginal dist of y | H1
  mu2 = pred_distr1$pred_mean
  var1 = pred_distr1$pred_var
  var2 = pred_distr1$pred_var
  Wass_dist = Wasserstein_distance_gp(mu1, mu2, var1, var2, dim = 1)
  return(Wass_dist)
}



###################
### one at time ###
###################

# # charge function evaluated at x
# q_gp = function(x_star, x_train, y_train, l01, type01){
#   pred_distr0 = getPredDistr(x_star, x_train, y_train, l01[1], type01[1])
#   pred_distr1 = getPredDistr(x_star, x_train, y_train, l01[2], type01[2])
#   mu1 = pred_distr0$pred_mean # mean of marginal dist of y | H1
#   mu2 = pred_distr1$pred_mean
#   var1 = pred_distr0$pred_cov
#   var2 = pred_distr1$pred_cov
#   Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2, dim = 1)
#   return(1.0 / Wass_dist)
# }
# 
# f_min_gp = function(candidate, D, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                  f0, f1, type, var_margy0, var_margy1, p, log_space = FALSE){
#   if(log_space == FALSE) {
#     result = q_gp(candidate, x_train, y_train, type01, l01)^k * sum(sapply(D, function(x_i) (q(x_i, x_train, y_train, type01, l01) / sqrt((x_i - candidate)^2))^k))
#     return(result)
#   } else{
#     # if has logSumExp library
#     terms = sapply(D, function(x_i) k * log(q(candidate, x_train, y_train, type01, l01)) + k * log(q(x_i, x_train, y_train, type01, l01)) - k * log(sqrt((x_i - candidate)^2)))
#     result = exp(logSumExp(terms))
#     return(result)
#   }
# }

# one at a time algorithm
MED_ms_gp_1pt = function(x_train, y_train, type01, l01,
                   N = 11, numCandidates = 10^5, k = 4, p = 2, 
                   xmin = 0, xmax = 1, nugget = NULL, genCandidates = 1){
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  # remove candidate points that are already in the design
  # this is necessary bc they lead to a posterior predictive distribution with
  # standard deviation/covariance = 0 and means equal to a design point
  # -> W = 0 -> division by 0 in q evaluation
  # candidates = candidates[!(candidates %in% x_train)]
  # but also can't have candidate points that are too close to them either...
  # lead to negative variances in predictive distribution somehow
  
  # -- Initialize 1st Design Point in D -- #
  w_evals = sapply(candidates, FUN = function(x) getW(x, x_train, y_train, type01, l01, nugget))
  xinitind = which.max(w_evals)
  return(candidates[xinitind])
}

