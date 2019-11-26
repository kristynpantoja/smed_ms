# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("variance_marginal_y.R")

##########
### 1D ###
##########

f_min = function(candidate, D, k, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                 f0, f1, type, var_margy0, var_margy1, p, alpha, buffer, log_space = FALSE){
  if(log_space == FALSE) {
    result = q(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
               f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^k * 
      sum(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                     f0, f1, type, var_margy0, var_margy1, p, alpha, buffer) / sqrt((x_i - candidate)^2))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                                f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)) + 
                     k * log(q(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, 
                               var_margy0, var_margy1, p, alpha, buffer)) - 
                     (k / 2) * log((x_i - candidate)^2))
    result = exp(logSumExp(terms))
    return(result)
  }
}

MED_ms_oneatatime = function(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                             f0 = NULL, f1 = NULL, type = NULL, N = 11, numCandidates = 10^5, k = 4, 
                             xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0, 
                             genCandidates = 1, initialpt = 1, var_margy0 = NULL, var_margy1 = NULL, 
                             jitter = FALSE, jittertype = 1, softmax = FALSE, threshold = FALSE, 
                             log_space = FALSE){
  # buffer = 0; genCandidates = 1; initialpt = 1; var_margy0 = NULL; var_margy1 = NULL 
  # jitter = FALSE; jittertype = 1; softmax = FALSE; threshold = FALSE; log_space = FALSE
  
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e - not necessary, if type is specified
  
  # some error checking, first:
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
  old_candidates = candidates
  # -- Initialize 1st Design Point in D -- #
  D = rep(NA, N)
  if(initialpt == 1){
    optimal_q = optimize(function(x) q(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, 
                                       type, var_margy0, var_margy1, p, alpha, buffer), interval = c(xmin, xmax))$minimum
    D[1] = optimal_q
  }
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    if(jitter == TRUE){
      candidates = old_candidates
      max_uniform = (xmax - xmin) / (numCandidates - 1)
      which_jitter_ind = 2:(numCandidates - 1)
      if(jittertype == 1) jitter_by = runif(length(which_jitter_ind), min = 0, max = max_uniform)
      if(jittertype == 2) jitter_by = runif(length(which_jitter_ind), min = 0, max = max_uniform / 2)
      candidates[which_jitter_ind] = candidates[which_jitter_ind] + jitter_by
    }
    # calculate cumulative TPE
    f_min_candidates = sapply(candidates, function(x) f_min(x, D[1:(i - 1)], k, mean_beta0, mean_beta1, 
                                                            var_beta0, var_beta1, var_e, f0, f1, 
                                                            type, var_margy0, var_margy1, p, alpha, buffer, log_space))
    f_opt = which.min(f_min_candidates)
    if(softmax == TRUE){
      # randomly select the new point
      # the probability is higher for candidates with lower TPEs
      sum_f_min_candidates = sum(f_min_candidates)
      softmax_candidates = 1 - (f_min_candidates / sum_f_min_candidates)
      D[i] = sample(candidates, 1, prob = softmax_candidates)
    } else if(threshold == TRUE){
      # randomly select from candidates with 5% of lowest cumulative TPEs
      order_ind_f_min_candidates = order(f_min_candidates)
      five_percent_smallest_ind = order(f_min_candidates)[1:(ceiling(numCandidates * 0.05))]
      xnew_ind = sample(five_percent_smallest_ind, 1)
      D[i] = candidates[xnew_ind]
    } else{
      D[i] = candidates[f_opt]
    }
  }
  return(D)
}





















##########
### 2D ###
##########

f_min_2d = function(candidate, D, k, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                    f0, f1, type, var_margy0, var_margy1, p, alpha, buffer, log_space = FALSE){
  if(log_space == FALSE) {
    result = q_2d(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                  f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^k * 
      sum(apply(D, 1, function(x_i) (q_2d(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                          f0, f1, type, var_margy0, var_margy1, p, alpha, buffer) / 
                                       sqrt(sum((x_i - candidate)^2)))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q_2d(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                                 f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)) + 
                     k * log(q_2d(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                  f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)) - 
                     k * log(sqrt((x_i - candidate)^2)))
    result = exp(logSumExp(terms))
    return(result)
  }
}

MED_ms_oneatatime_2d = function(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                     f0 = NULL, f1 = NULL, type = NULL, N = 11, numCandidates = 10^5, k = 4, 
                     xmin = 0, xmax = 1, p = 3, alpha = NULL, buffer = 0, 
                     log_space = FALSE, genCandidates = 1, var_margy0 = NULL, var_margy1 = NULL){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
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
  w_evals = apply(candidates, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                    var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
                                                                    var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
  xinitind = which.max(w_evals)
  D[1, ] = candidates[xinitind, ]
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, D[1:(i - 1), , drop = FALSE], k, mean_beta0, mean_beta1, 
                                                                 var_beta0, var_beta1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
    #print(i)
  }
  return(D)
}


