# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("variance_marginal_y.R")

##########
### 1D ###
##########

f_min_nodata = function(candidate, D, k, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                        f0, f1, type, var_margy0, var_margy1, p, alpha, buffer, log_space = FALSE){
  if(log_space == FALSE) {
    result = q_nodata(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
               f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^k * 
      sum(sapply(D, function(x_i) (q_nodata(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
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

generate_MMED_nodata_oneatatime = function(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
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
    optimal_q = optimize(function(x) q_nodata(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, 
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
    f_min_candidates = sapply(candidates, function(x) f_min_nodata(x, D[1:(i - 1)], k, mean_beta0, mean_beta1, 
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





###


f_min_fast_nodata = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                             f0, f1, type, var_margy0, var_margy1, p, alpha, buffer){
  
  q(candidate_jk, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
    f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^gamma_k * 
    max(sapply(D_k, function(x_i) q(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                    f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^gamma_k / 
                 sqrt((x_i - candidate_jk)^2)))
}

generate_MMED_nodata_fast = function(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                     f0 = NULL, f1 = NULL, type = NULL, N = 11, K = 10, 
                                     xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0, numCandidates = NULL, 
                                     genCandidates = 1, initialpt = 1, var_margy0 = NULL, var_margy1 = NULL){
  
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  if(is.null(numCandidates)) numCandidates = N
  
  # -- Create hypothesized models -- #
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
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = rep(NA, numCandidates)
  if(genCandidates == 1) C1 = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) C1 = sort(runif(numCandidates, xmin, xmax))
  if(genCandidates == 3) C1 = mined::Lattice(numCandidates, p = p)
  D1 = C1
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K = 1, return the space-filling design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = N, ncol = K)
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = c(1:(K - 1)) / (K - 1) # Shouldn't this be 0:(K-1)/(k-1) ? ************************************
  # -- the final step should be gamma = 1 because then we optimize the correct criterion
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  optimal_q = optimize(function(x) q_nodata(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, 
                                            type, var_margy0, var_margy1, p, alpha, buffer), interval = c(xmin, xmax))$minimum
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    
    if(initialpt == 2){
      # Get candidates in neighborhood L1k = (lower, upper):
      # for j = 1
      # get candidates in L_1k
      R1k = min(abs(D[-1, k] - D[1, k])) # radius of L1k
      lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ************************************
      upper = min(D[1, k] + R1k, 1) # HERE IT IS BECAUSE o/w GET NaNs in q evaluation! why, though/? *******************
      # candidates from space-filling design, tildeD1_kplus1
      # In general the number of local points does no need to be N ************************************ans
      # so I suggest introducing a N_L. You can set N_L = N for now ************************************ans
      # but we may decide to change it later. ************************************ans
      tildeD1_kplus = rep(NA, numCandidates)
      if(genCandidates == 1) tildeD1_kplus = seq(from = lower, to = upper, length.out = numCandidates)
      if(genCandidates == 2) tildeD1_kplus = runif(numCandidates, lower, upper)
      if(genCandidates == 3) tildeD1_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
      # save the candidates to be used in future designs
      C[[1]] = c(C[[1]], tildeD1_kplus)
      # criterion to choose first candidate from candidate set: 
      # the point at which f1 and f2 are most different
      w_evals = sapply(C[[1]], FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                      var_marginaly(x, var_e, var_beta0, type, var_margy0), 
                                                                      var_marginaly(x, var_e, var_beta1, type, var_margy1)))
      # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
      xinitind = which.max(w_evals)
      
      D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
      
      #Wass_D[1, k] = ... # a next step would be to matricize/preallocate these values for faster computing! ********************
      
    } else{
      D[1, k + 1] = optimal_q
    }
    
    
    # for j = 2:N
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
      lower = max(D[j, k] - R_jk, 0) 
      upper = min(D[j, k] + R_jk, 1)
      tildeDj_kplus = rep(NA, numCandidates)
      if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
      if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
      if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
      C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
      
      f_min_candidates = sapply(C[[j]], function(x) f_min_nodata_fast(x, D[1:(j - 1), k + 1], gammas[k], 
                                                                      mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                                                                      f0, f1, type, var_margy0, var_margy1, p, alpha, buffer))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k + 1] = C[[j]][chosen_cand]
      
    }
  }
  
  return(list("D" = D, "candidates" = C))
}




















##########
### 2D ###
##########

f_min_nodata_2d = function(candidate, D, k, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                    f0, f1, type, var_margy0, var_margy1, p, alpha, buffer, log_space = FALSE){
  if(log_space == FALSE) {
    result = q_nodata_2d(candidate, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                  f0, f1, type, var_margy0, var_margy1, p, alpha, buffer)^k * 
      sum(apply(D, 1, function(x_i) (q_nodata_2d(x_i, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
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

MMED_nodata_oneatatime_2d = function(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
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
    f_min_candidates = apply(candidates, 1, function(x) f_min_nodata_2d(x, D[1:(i - 1), , drop = FALSE], k, mean_beta0, mean_beta1, 
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


