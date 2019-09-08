require("wasserstein_distance.R")
require("charge_function_q.R")
require("variance_marginal_y.R")

##########
### 1D ###
##########

f_min_fast = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                      f0, f1, type, var_margy0, var_margy1, p){
  
  q(candidate_jk, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
    f0, f1, type, var_margy0, var_margy1, p)^gamma_k * 
    max(sapply(D_k, function(x_i) q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                    f0, f1, type, var_margy0, var_margy1, p)^gamma_k / sqrt((x_i - candidate_jk)^2)))
}

MED_ms_fast = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                       N = 11, numCandidates = NULL, K = 10, p = 2, xmin = 0, xmax = 1, 
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
      R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
      lower = max(D[j, k] - R_jk, 0) 
      upper = min(D[j, k] + R_jk, 1)
      tildeDj_kplus = rep(NA, numCandidates)
      if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
      if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
      if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
      C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
      
      f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], 
                                                               mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                               f0, f1, type, var_margy0, var_margy1, p))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k + 1] = C[[j]][chosen_cand]
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}



##########
### 2D ###
##########

f_min_fast_2d = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                         f0, f1, type, var_margy0, var_margy1, p){
  q_2d(candidate_jk, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
       f0, f1, type, var_margy0, var_margy1, p)^gamma_k * 
    max(apply(D_k, 1, function(x_i) q_2d(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                         f0, f1, type, var_margy0, var_margy1, p)^gamma_k / sqrt(sum((x_i - candidate_jk)^2))))
}

MED_ms_fast_2d = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                          f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                          N = 11, numCandidates = NULL, K = 10, p = 2, xmin = 0, xmax = 1, seed = 1){
  
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
  
  # -- Make D_1, space-filling design -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  sqrtNumCand = NULL
  if(is.null(numCandidates)){
    # make the number of candidates approximately equal to N (must be an integer square root, to make a grid)
    sqrtNumCand = floor(sqrt(N))
  }
  if(length(numCandidates) == 1){ # 
    sqrtNumCand = floor(sqrt(numCandidates)) # to be able to make a grid of candidates
  }
  numCandidates = c(sqrtNumCand, sqrtNumCand)
  numCandidatesTtl = sqrtNumCand^2
  candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
  candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
  C1 = cbind(rep(candidates_x1, each = numCandidates[2]), 
             rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
  spots_left = N - dim(C1)[1]
  if(spots_left != 0){
    set.seed(seed)
    candidates_x1_leftover = runif(spots_left, xmin, xmax)
    candidates_x2_leftover = runif(spots_left, xmin, xmax)
    D1 =  rbind(C1, cbind(candidates_x1_leftover, candidates_x2_leftover))
  } else{
    D1 = C1
  }
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  # if do this instead of calculating q for each pair each time, may save time and space
  
  # -- If K = 1, return the space-filling design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  Dkplus = matrix(rep(NA, N), N, 2)
  Dk = D1
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = c(1:(K - 1)) / (K - 1)
  # -- the final step should be gamma = 1 because then we optimize the correct criterion
  # save candidates for each K
  C <- array(NA, dim = c(numCandidatesTtl * K, 2, N))
  for (j in 1:N){
    C[ , , j] = rbind(C1, matrix(rep(NA, 2 * numCandidatesTtl * (K - 1)), numCandidatesTtl * (K - 1), 2))
  }
  
  # at index k, determine the next design k + 1
  
  ## For j = 1, i.e. 1st point in design k + 1: it's always the same point for all k ##
  # criterion to choose first candidate from candidate set: 
  # the point at which f1 and f2 are most different
  w_evals = apply(C[,,1], 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                    var_marginaly_2d(as.vector(x), var_mean0, var_e, type[1], var_margy0), 
                                                                    var_marginaly_2d(as.vector(x), var_mean1, var_e, type[2], var_margy1)))
  # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
  xinitind = which.max(w_evals)
  
  for(k in 1:(K - 1)){
    Dkplus[1, ] = C[,,1][xinitind, ] # x1, first element of set of design points, D
    
    # for j = 2:N
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      diff = (Dk[-j, ] - matrix(rep(Dk[j, ], N - 1), N - 1, 2, byrow = TRUE))^2
      Rjk = min(sqrt(rowSums(diff)))
      lower_x1 = max(Dk[j, 1] - Rjk, 0); lower_x2 = max(Dk[j, 2] - Rjk, 0)
      upper_x1 = min(Dk[j, 1] + Rjk, 1); upper_x2 = min(Dk[j, 2] + Rjk, 1)
      tildeDj_kplus_x1 = seq(from = lower_x1, to = upper_x1, length.out = numCandidates[1])
      tildeDj_kplus_x2 = seq(from = lower_x2, to = upper_x2, length.out = numCandidates[2])
      tildeDj_kplus = cbind(rep(tildeDj_kplus_x1, each = numCandidates[2]), 
                            rep(tildeDj_kplus_x2, times = numCandidates[1]))
      C[(k * numCandidatesTtl + 1):((k + 1) * numCandidatesTtl) , , j] = tildeDj_kplus # This is now C_j^{k+1}
      
      f_min_candidates = apply(C[ , , j], 1, function(x) f_min_fast_2d(x, Dkplus[1:(j - 1), , drop = FALSE], gammas[k], 
                                                                       mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                       f0, f1, type, var_margy0, var_margy1, p))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      Dkplus[j, ] = C[ , , j][chosen_cand, ]
      #print(paste("k:", k, " out of ", (K - 1)," --- j:", j, " out of ", N))
    }
    Dk = Dkplus
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = Dkplus, "candidates" = C))
}
