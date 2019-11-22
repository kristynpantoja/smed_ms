# require("wasserstein_distance.R")
# require("charge_function_q.R")

# 1d Covariance functions
# radial aka squared exponential aka gaussian
# phi_sqexp = function(Xi,Xj, l) exp(-0.5 * ((Xi - Xj) / ls)^2)
phi_sqexp = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
# exponential, a.k.a. matern, nu = 1/2
phi_exp = function(Xi, Xj, l) exp(- abs(Xi - Xj) / l)
# matern, nu = 3/2
phi_matern2 = function(Xi, Xj, l){
  r = abs(Xi - Xj)
  (1 + (sqrt(3) * r / l)) * exp(- sqrt(3) * r / l)
}
# general matern, l = 1
phi_matern = function(Xi, Xj, l, nu = 3/2) RandomFieldsUtils::matern(abs(Xi - Xj), nu, scaling = "matern")
# periodic, but a generic version
# phi_periodic = function(Xi, Xj, l){
#   r = Xi - Xj
#   exp(-2 * (sin(r / 2) / l)^2)
# }
# a periodic kernel that allows to adjust the period, p. 
phi_periodic = function(Xi, Xj, l, p = pi / 24){
  r = abs(Xi - Xj)
  sinsq_arg = pi * r / p
  exp(-2 * (sin(sinsq_arg / 2) / l)^2)
}

# compute covariance matrix
getCov = function(X, Y, type = 1, l = 1){ # change l = 1 to a list theta = NULL
  phi = NULL
  if(type == 1) phi = phi_sqexp
  else if(type == 2) phi = phi_matern
  else if(type == 3) phi = phi_exp
  else if(type == 4) phi = phi_matern2
  else if(type == 5) phi = phi_periodic
  else stop("invalid type specification of covariance function for GP")
  outer(X, Y, FUN = phi, l)
}

getPredDistr = function(x_star, x_train, y_train, type, l, nugget = NULL){
  k_star = t(as.matrix(getCov(x_star, x_train, type, l)))
  if(is.null(nugget)) K_obs = getCov(x_train, x_train, type, l)
  else K_obs = getCov(x_train, x_train, type, l) + diag(rep(nugget, length(x_train)))
  pred_mean = t(k_star) %*% solve(K_obs, y_train)
  pred_cov = getCov(x_star, x_star, type, l) - t(k_star) %*% solve(K_obs, k_star)
  return(list("pred_mean" = as.vector(pred_mean), "pred_var" = pred_cov))
}


getPredDistrSeq = function(x_seq, x_train, y_train, type, l, nugget = NULL){
  k_star = t(getCov(x_seq, x_train, type, l))
  if(is.null(nugget)) K_obs = getCov(x_train, x_train, type, l)
  else K_obs = getCov(x_train, x_train, type, l) + diag(rep(nugget, length(x_train)))
  pred_mean = t(k_star) %*% solve(K_obs, y_train)
  pred_cov = getCov(x_seq, x_seq, type, l) - (t(k_star) %*% solve(K_obs, k_star))
  return(list("pred_mean" = as.vector(pred_mean), "pred_var" = pred_cov))
}

getW = function(x_star, x_train, y_train, type01, l01, nugget){
  pred_distr0 = getPredDistr(x_star, x_train, y_train, type01[1], l01[1], nugget)
  pred_distr1 = getPredDistr(x_star, x_train, y_train, type01[2], l01[2], nugget)
  mu1 = pred_distr0$pred_mean # mean of marginal dist of y | H1
  mu2 = pred_distr1$pred_mean
  var1 = pred_distr0$pred_var
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





###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# using data (posterior predictive distribution of y)

f_min_data_gp = function(candidate, D, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer){
  result = q_data_gp(candidate, Kinv0, Kinv1, initD, y, var_e, type, l, p, 
                     alpha, buffer)^k * 
    sum(sapply(D, function(x_i) (q_data_gp(x_i, Kinv0, Kinv1, initD, y, var_e, type, l, p, 
                                           alpha, buffer) / sqrt((x_i - candidate)^2))^k))
  return(result)
}


add_MED_ms_oneatatime_data_gp_old = function(initD, y, type, l, var_e, N2 = 11, numCandidates = 10^5, k = 4, p = 2, 
                                             xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                             genCandidates = 1, candidates = NULL){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta relies on Kinv for each hypothesized kernel function
  if(is.null(nugget)){
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  } else{
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, length(initD))))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, length(initD))))
  }
  
  # check for any bad initial points in initD
  w_initD = sapply(initD, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  if(length(which(w_initD == 0)) != 0){
    initD = initD[-which(w_initD == 0)]
    y = y[-which(w_initD == 0)]
    w_initD = w_initD[-which(w_initD == 0)]
    if(is.null(nugget)){
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, length(initD))))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, length(initD))))
    }
  }
  
  initN = length(initD)
  ttlN = initN + N2
  
  
  # -- Generate Candidate Points -- #
  if(!is.null(candidates)){
    numCandidates = length(candidates)
  } else{
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  # remove candidate points that are already in the design
  # this is necessary bc they lead to a posterior predictive distribution with
  # standard deviation/covariance = 0 and means equal to a design point
  # -> W = 0 -> division by 0 in q evaluation
  # candidates = candidates[!(candidates %in% x_train)] # actually not necessray, 
  # considering our edit Wasserstein_distance_gp = 0 when can't take sqrt
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  # w_evals = sapply(candidates, FUN = function(x) getW(x, initD, y, type, l, nugget))
  w_evals = sapply(candidates, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex]
  is_xopt_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_xopt_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data_gp(x, initD, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = xoptindex
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data_gp(x, c(initD, D[1:(i - 1)]), Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
    D_ind[i] = f_opt
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD, 
              "candidates" = candidates, "indices" = D_ind))
}

add_MED_ms_oneatatime_data_gp = function(initD, y, type, l, var_e, N2 = 11, numCandidates = 10^5, k = 4, p = 2, 
                                         xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                         genCandidates = 1, candidates = NULL){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta
  if(is.null(nugget)){
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  } else{
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, length(initD))))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, length(initD))))
  }
  
  # check for any bad initial points in initD
  # w_initD = sapply(initD, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  # if(length(which(w_initD == 0)) != 0){
  #   initD = initD[-which(w_initD == 0)]
  #   y = y[-which(w_initD == 0)]
  #   w_initD = w_initD[-which(w_initD == 0)]
  #   if(is.null(nugget)){
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  #   } else{
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, length(initD))))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, length(initD))))
  #   }
  # }
  
  initN = length(initD)
  ttlN = initN + N2
  
  
  # -- Generate Candidate Points -- #
  if(!is.null(candidates)){
    numCandidates = length(candidates)
  } else{
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  # remove candidate points that are already in the design
  # this is necessary bc they lead to a posterior predictive distribution with
  # standard deviation/covariance = 0 and means equal to a design point
  # -> W = 0 -> division by 0 in q evaluation
  # candidates = candidates[!(candidates %in% x_train)] # actually not necessray, 
  # considering our edit Wasserstein_distance_gp = 0 when can't take sqrt
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  w_evals = sapply(candidates, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex]
  is_xopt_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_xopt_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data_gp(x, initD, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = xoptindex
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data_gp(x, D[1:(i - 1)], Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
    D_ind[i] = f_opt
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD, 
              "candidates" = candidates, "indices" = D_ind))
}





