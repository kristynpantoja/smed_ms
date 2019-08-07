


################################################
### one-at-a-time method with fast criterion ###
################################################

# this is just fast algorithm with k = 1, but using max instead of sum over the design points

f_min_mix1 = function(candidate, D, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  q(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p) * 
    max(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p) / sqrt((x_i - candidate)^2))))
}

# sequentially pick each point (using a greedy algorithm) but by minimizing the fast
#  algorithm's criterion instead (max instead of summation)
MED_ms_mix1 = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                       N = 11, numCandidates = 10^5, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
                       genCandidates = 1, initialpt = 1){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
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
  
  # -- Initialize 1st Design Point in D -- #
  D = rep(NA, N)
  if(initialpt == 2){
    xinitind = which.max(abs(f0(candidates) - f1(candidates)))
    D[1] = candidates[xinitind]
  } else{
    D[1] = optimize(function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, 
                                  type, var_margy0, var_margy1, p), 
                    interval = c(xmin, xmax))$minimum
  }
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_mix1(x, D[1:(i - 1)], mean_beta0, mean_beta1, 
                                                                 var_mean0, var_mean1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
  }
  return(D)
}











################################################
### fast method with one-at-a-time criterion ###
################################################

# have points converge to those which minimize the one-at-a-time algorithm's 
#  criterion (summation raised to k = 4 power) over K stages

f_min_mix2 = function(candidate_jk, D_k, gamma_k, k_power, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                      f0, f1, type, var_margy0, var_margy1, p){
  q(candidate_jk, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
    f0, f1, type, var_margy0, var_margy1, p)^(gamma_k * k_power) * 
    sum(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                     f0, f1, type, var_margy0, var_margy1, p)^gamma_k / sqrt((x_i - candidate_jk)^2))^k_power))
}

MED_ms_mix2 = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                       N = 11, numCandidates = NULL, K = 10, k_power = 4, p = 2, xmin = 0, xmax = 1, 
                       genCandidates = 1, initialpt = 1){
  
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
  
  optimal_q = optimize(function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, 
                                     type, var_margy0, var_margy1, p), interval = c(xmin, xmax))$minimum
  
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
                                                                      var_marginaly(x, var_e, var_mean0, type, var_margy0), 
                                                                      var_marginaly(x, var_e, var_mean1, type, var_margy1)))
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
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # Shouldn't this be -R_jk ... and 0 not 1? ************************************fixed
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        # Did you check their code to verify this is how to construct the candidate set? ************************************
        # (the paper doesn't seem very clear) ************************************
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_mix2(x, D[1:(j - 1), k + 1], gammas[k], k_power, 
                                                                 mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                 f0, f1, type, var_margy0, var_margy1, p))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_mix2(x, D[1:(j - 1), k + 1], gammas[k], k_power, 
                                                                 mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                 f0, f1, type, var_margy0, var_margy1, p))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}