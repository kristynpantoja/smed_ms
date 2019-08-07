

#####################
#####################
#####################
#####################
###### 2-DIM'L ######
#####################
#####################
#####################
#####################

# Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly_2d = function(x, var_mean, var_e, type, var_margy){
  # type:
  #   1 for linear model without slope
  #   2 for linear model with slope
  #   3 for quadratic model with slope
  #   4 for two-dimentional model
  if(!is.null(type)){
    if(type == 4) var_mean[1] + x[1]^2 * var_mean[2] + x[2]^2 * var_mean[3] + var_e
    else stop(paste("invalid type given : ", type))
  } else{
    if(!is.null(var_margy)) var_margy(x = x, var_mean = var_mean, var_e = var_e)
    else stop("cannot compute var_marginaly: no type given, and no var_margy fn given")
  }
}

# charge function evaluated at x
q_2d = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly_2d(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly_2d(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(1.0 / Wass_dist^(1 / (2 * p)))
}


###################
### one at time ###
###################

f_min_2d = function(candidate, D, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                    f0, f1, type, var_margy0, var_margy1, p, log_space = FALSE){
  if(log_space == FALSE) {
    
    result = q_2d(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)^k * sum(apply(D, 1, function(x_i) (q_2d(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p) / sqrt(sum((x_i - candidate)^2)))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q_2d(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) + k * log(q_2d(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) - k * log(sqrt((x_i - candidate)^2)))
    result = exp(logSumExp(terms))
    return(result)
  }
}

MED_ms_2d = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                     f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                     N = 11, numCandidates = 10^5, k = 4, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
                     genCandidates = 1){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
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
  f0_cand = apply(candidates, 1, FUN = f0)
  f1_cand = apply(candidates, 1, FUN = f1)
  xinitind = which.max(abs(f0_cand - f1_cand)) # same index as max wasserstein distance, but maybe at some point implement that one instead
  D[1, ] = candidates[xinitind, ]
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, D[1:(i - 1), , drop = FALSE], k, mean_beta0, mean_beta1, 
                                                                 var_mean0, var_mean1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
  }
  return(D)
}


############
### fast ###
############


f_min_fast_2d = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                         f0, f1, type, var_margy0, var_margy1, p){
  q_2d(candidate_jk, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
       f0, f1, type, var_margy0, var_margy1, p)^gamma_k * 
    max(apply(D_k, 1, function(x_i) q_2d(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                         f0, f1, type, var_margy0, var_margy1, p)^gamma_k / sqrt(sum((x_i - candidate_jk)^2))))
}

MED_ms_fast_2d = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                          f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                          N = 11, numCandidates = NULL, K = 10, p = 2, xmin = 0, xmax = 1){
  
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  if(is.null(numCandidates)) numCandidates = N
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Make D_1, space-filling design -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  if(length(numCandidates) == 1){
    sqrtN = floor(sqrt(numCandidates))
    numCandidates = c(sqrtN, sqrtN)
    N = sqrtN^2
  }
  candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
  candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
  C1 = cbind(rep(candidates_x1, each = numCandidates[2]), 
             rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
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
  Dkplus = matrix(rep(NA, N), N, 2)
  Dk = D1
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = c(1:(K - 1)) / (K - 1)
  # -- the final step should be gamma = 1 because then we optimize the correct criterion
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1: ##
    # Get candidates in neighborhood L1k = (lower, upper):
    # for j = 1
    # get candidates in L_1k
    diff = (Dk[-1, ] - matrix(rep(Dk[1, ], N - 1), N - 1, 2, byrow = TRUE))^2
    R1k = min(sqrt(rowSums(diff))) # minimum euclidean distance b/t
    lower = max(Dk[1, ] - R1k, 0)
    upper = min(Dk[1, ] + R1k, 1)
    tildeD1_kplus_x1 = seq(from = lower, to = upper, length.out = numCandidates[1])
    tildeD1_kplus_x2 = seq(from = lower, to = upper, length.out = numCandidates[2])
    tildeD1_kplus = cbind(rep(tildeD1_kplus_x1, each = numCandidates[2]), 
                          rep(tildeD1_kplus_x2, times = numCandidates[1]))
    # save the candidates to be used in future designs
    C[[1]] = rbind(C[[1]], tildeD1_kplus)
    # criterion to choose first candidate from candidate set: 
    # the point at which f1 and f2 are most different
    w_evals = apply(C[[1]], 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                      var_marginaly_2d(as.vector(x), var_mean0, var_e, type[1], var_margy0), 
                                                                      var_marginaly_2d(as.vector(x), var_mean1, var_e, type[2], var_margy1)))
    # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
    xinitind = which.max(w_evals)
    Dkplus[1, ] = C[[1]][xinitind, ] # x1, first element of set of design points, D
    
    # for j = 2:N
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      diff = (Dk[-j, ] - matrix(rep(Dk[j, ], N - 1), N - 1, 2, byrow = TRUE))^2
      Rjk = min(sqrt(rowSums(diff)))
      lower = max(Dk[j, ] - Rjk, 0) 
      upper = min(Dk[j, ] + Rjk, 1)
      tildeDj_kplus_x1 = seq(from = lower, to = upper, length.out = numCandidates[1])
      tildeDj_kplus_x2 = seq(from = lower, to = upper, length.out = numCandidates[2])
      tildeDj_kplus = cbind(rep(tildeDj_kplus_x1, each = numCandidates[2]), 
                            rep(tildeDj_kplus_x2, times = numCandidates[1]))
      C[[j]] = rbind(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
      
      f_min_candidates = apply(C[[j]], 1, function(x) f_min_fast_2d(x, Dkplus[1:(j - 1), , drop = FALSE], gammas[k], 
                                                                    mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                    f0, f1, type, var_margy0, var_margy1, p))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      Dkplus[j, ] = C[[j]][chosen_cand, ]
    }
    Dk = Dkplus
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = Dkplus, "candidates" = C))
}


##################
### evaluate D ###
##################

### --- Compute Criteria --- ###

totalPE_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2))
      counter = counter + 1
    }
  }
  return(sum(pairwise_PEs))
}


crit_1atatime_2d = function(D, N, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  if(log_space == FALSE) {
    numPairs = N * (N - 1) / 2
    pairwise_PEs = rep(NA, numPairs)
    counter = 1
    qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_PEs[counter] = (qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2)))^k
        counter = counter + 1
      }
    }
    return((sum(pairwise_PEs))^(1/k))
  } else{
    if(N != dim(D)[1]) stop("N is not the same as length of D")
    numPairs = N * (N - 1) / 2
    pairwise_terms = rep(NA, numPairs)
    counter = 1
    logqD = apply(D, 1, FUN = function(x) log(q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                   f0, f1, type, var_margy0, var_margy1, p)))
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_terms[counter] = k * logqD[i] + k * logqD[j] - k * log(sqrt(sum((D[i, ] - D[j, ])^2)))
        counter = counter + 1
      }
    }
    return(exp((1 / k) * logSumExp(pairwise_terms)))
  }
}

crit_fast_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2))
      counter = counter + 1
    }
  }
  return(max(pairwise_PEs))
}



