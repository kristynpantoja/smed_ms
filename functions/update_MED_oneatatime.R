# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("variance_marginal_y.R")
# require("construct_design_matrix.R")
# require("posterior_variance.R")
# require("posterior_mean.R")

##########
### 1D ###
##########

# without using data (marginal distribution of y)

add_MED_ms_oneatatime = function(initD, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                     f0 = NULL, f1 = NULL, type = NULL, N2 = 11, numCandidates = 10^5, k = 4, 
                                     xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0, 
                                     wasserstein0 = 1, genCandidates = 1, 
                                     var_margy0 = NULL, var_margy1 = NULL, log_space = FALSE){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  if(wasserstein0 == 1){
    w_initD = sapply(initD, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                     var_marginaly(x, var_beta0, var_e, type[1], var_margy0), 
                                                                     var_marginaly(x, var_beta1, var_e, type[2], var_margy1)))
    if(length(which(w_initD == 0)) != 0){
      initD = initD[-which(w_initD == 0)]
      w_initD = w_initD[-which(w_initD == 0)]
    }
  }
  if(wasserstein0 == 2){
    if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 0.0001")
    buffer == 0.0001
  }
  
  # other variables and checks
  initN = length(initD)
  ttlN = initN + N2
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 1) f0 = function(x) mean_beta0 * x
    else if(type[1] == 2) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x
    else if(type[1] == 3) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x + mean_beta0[3] * x^2
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[2] == 1) f1 = function(x) mean_beta1 * x
    else if(type[2] == 2) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x
    else if(type[2] == 3) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  # wass_dist_fun = function(x){
  #   distr1.mu = f0(x)
  #   distr1.var = var_marginaly(x, var_beta0, var_e, type = type[1], var_margy0)
  #   distr2.mu = f1(x)
  #   distr2.var = var_marginaly(x, var_beta1, var_e, type = type[2], var_margy1)
  #   wass_dist = Wasserstein_distance(distr1.mu, distr2.mu, distr1.var, distr2.var)
  #   return(wass_dist)
  # }
  # max.sep.loc2 = optimize(wass_dist_fun, lower=xmin, upper=xmax, maximum = TRUE)$maximum
  # D[1] = max.sep.loc2
  # this does the same thing: but may want to do the above if we want to save wasserstein distances
  # to make the code run raster at some point. for now, do this.
  optimal_q = optimize(function(x) q(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, 
                                     type, var_margy0, var_margy1, p, alpha, buffer), interval = c(xmin, xmax))$minimum
  xopt = optimal_q
  is_x_max_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min(x, initD, k, mean_beta0, mean_beta1, 
                                                            var_beta0, var_beta1, var_e, f0, f1, 
                                                            type, var_margy0, var_margy1, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
  } else{
    D[1] = xopt
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min(x, c(initD, D[1:(i - 1)]), k, mean_beta0, mean_beta1, 
                                                            var_beta0, var_beta1, var_e, f0, f1, 
                                                            type, var_margy0, var_margy1, p, alpha, buffer, log_space))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
  }
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD))
}


# using data (posterior predictive distribution of y)

f_min_data2 = function(candidate, D, postmean0, postmean1, postvar0, postvar1, var_e, type, p, alpha, buffer){
  result = q_data2(candidate, postmean0, postmean1, postvar0, postvar1, var_e, type, p, 
                  alpha, buffer)^k * 
    sum(sapply(D, function(x_i) (q_data2(x_i, postmean0, postmean1, postvar0, postvar1, var_e, type, p, 
                                        alpha, buffer) / sqrt((x_i - candidate)^2))^k))
  return(result)
}

add_MED_ms_oneatatime_data = function(initD, y, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                      f0 = NULL, f1 = NULL, type = NULL, N2 = 11, numCandidates = 10^5, k = 4, 
                                      xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0, 
                                      wasserstein0 = 1, genCandidates = 1, 
                                      var_margy0 = NULL, var_margy1 = NULL, log_space = FALSE){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta
  postvar0 = postvar(initD, initN, var_e, var_beta0, type[1])
  postmean0 = postmean(y, initD, initN, mean_beta0, var_beta0, var_e, type[1])
  postvar1 = postvar(initD, initN, var_e, var_beta1, type[2])
  postmean1 = postmean(y, initD, initN, mean_beta1, var_beta1, var_e, type[2])
  
  if(wasserstein0 == 1){
    w_initD = sapply(initD, FUN = function(x) Wasserstein_distance_postpred2(x, postmean0, postmean1, postvar0, postvar1, var_e, type))
    if(length(which(w_initD == 0)) != 0){
      initD = initD[-which(w_initD == 0)]
      w_initD = w_initD[-which(w_initD == 0)]
      y = y[-which(w_initD == 0)]
      postvar0 = postvar(initD, initN, var_e, var_beta0, type[1])
      postmean0 = postmean(y, initD, initN, mean_beta0, var_beta0, var_e, type[1])
      postvar1 = postvar(initD, initN, var_e, var_beta1, type[2])
      postmean1 = postmean(y, initD, initN, mean_beta1, var_beta1, var_e, type[2])
    }
  }
  if(wasserstein0 == 2){
    if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 0.0001")
    buffer == 0.0001
  }
  
  # other variables and checks
  initN = length(initD)
  ttlN = initN + N2
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 1) f0 = function(x) mean_beta0 * x
    else if(type[1] == 2) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x
    else if(type[1] == 3) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x + mean_beta0[3] * x^2
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[2] == 1) f1 = function(x) mean_beta1 * x
    else if(type[2] == 2) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x
    else if(type[2] == 3) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  # wass_dist_fun = function(x){
  #   distr1.mu = f0(x)
  #   distr1.var = var_marginaly(x, var_beta0, var_e, type = type[1], var_margy0)
  #   distr2.mu = f1(x)
  #   distr2.var = var_marginaly(x, var_beta1, var_e, type = type[2], var_margy1)
  #   wass_dist = Wasserstein_distance(distr1.mu, distr2.mu, distr1.var, distr2.var)
  #   return(wass_dist)
  # }
  # max.sep.loc2 = optimize(wass_dist_fun, lower=xmin, upper=xmax, maximum = TRUE)$maximum
  # D[1] = max.sep.loc2
  # this does the same thing: but may want to do the above if we want to save wasserstein distances
  # to make the code run raster at some point. for now, do this.
  optimal_q = optimize(function(x) q_data2(x, postmean0, postmean1, postvar0, postvar1, var_e, type, p, 
                                           alpha, buffer), interval = c(xmin, xmax))$minimum
  xopt = optimal_q     
  is_x_max_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data2(x, initD, postmean0, postmean1, postvar0, postvar1, var_e, type, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
  } else{
    D[1] = xopt
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_data2(x, c(initD, D[1:(i - 1)]), postmean0, postmean1, postvar0, postvar1, var_e, type, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD))
}




































##########
### 2D ###
##########

# hasn't been tested
add_MED_ms_oneatatime_2d = function(initD, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                        f0 = NULL, f1 = NULL, type = NULL, N = 11, numCandidates = 10^5, k = 4, 
                                        xmin = 0, xmax = 1, p = 3, alpha = NULL, buffer = 0, y = NULL, 
                                        log_space = FALSE, genCandidates = 1, wasserstein0 = 1, 
                                        var_margy0 = NULL, var_margy1 = NULL){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  if(wasserstein0 == 1){
    w_initD = apply(initD, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                     var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
                                                                     var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
    initD = initD[-which(w_initD == 0), ]
    w_initD = w_initD[-which(w_initD == 0), ]
  }
  if(wasserstein0 == 2){
    if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 0.0001")
    buffer == 0.0001
  }
  # other variables and checks
  initN = dim(initD)[1]
  ttlN = initN + N
  if(dim(initD)[2] != 2) stop("initD is not a matrix of size initN x 2")
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
    else if(type[1] == 5) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[1]^2 + mean_beta0[4] * x[2] + mean_beta0[3] * x[2]^2
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
    else if(type[1] == 5) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[1]^2 + mean_beta1[4] * x[2] + mean_beta1[3] * x[2]^2
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Generate Candidate Points -- #
  
  if(length(numCandidates) == 1) numCandidates = c(floor(sqrt(numCandidates)), floor(sqrt(numCandidates)))
  # sometimes one factor is more important, so gets more candidates, than the other factor.
  if(genCandidates == 1){
    candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
    candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
    candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
                       rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
  }
  if(genCandidates == 2){
    candidates_x1 = sort(runif(numCandidates[1], min = xmin, max = xmax))
    candidates_x2 = sort(runif(numCandidates[2], min = xmin, max = xmax))
    candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
                       rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
  }
  
  # -- Initialize 1st Design Point in D -- #
  D = matrix(rep(NA, N * 2), N, 2)
  w_cand = apply(candidates, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                       var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
                                                                       var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
  xinitind = which.max(w_cand)
  xmaxW = candidates[xinitind, ]
  is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), xmaxW)) # give tolerance?
  if(is_x_max_in_initD){
    N = N - 1
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, initD, k, mean_beta0, mean_beta1, 
                                                                 var_beta0, var_beta1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[1, ] = xnew
  } else{
    D[1, ] = x_maxW
  }
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, rbind(initD, D[1:(i - 1), , drop = FALSE]), k, mean_beta0, mean_beta1, 
                                                                 var_beta0, var_beta1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
    print(i)
  }
  return(list("initD" = old_initD, "addD" = D, "updatedD" = rbind(old_initD, D), "q_initD" = initD))
}

