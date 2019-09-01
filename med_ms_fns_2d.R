


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
    else if(type == 5) var_mean[1] + x[1]^2 * var_mean[2] + x[1]^4 * var_mean[3] + x[2]^2 * var_mean[4] + x[2]^4 * var_mean[5] + var_e
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
  return(1.0 / Wass_dist)
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
    #print(i)
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



##################################################################
##################################################################
# ADDING POINTS TO A DESIGN #
##################################################################
##################################################################

# IGNORING POINTS THAT GIVE WASSERSTEIN DISTANCE OF 0!!!! IS THIS THE RIGHT THING TO DO?
# could instead add a buffer to the wasserstein distance? mu0 - mu1 + 0.01 * (xmax - xmin) or something?
MED_ms_2d_add = function(initD, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                         f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                         N = 11, numCandidates = 10^5, k = 4, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
                         genCandidates = 1){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  w_initD = apply(initD, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                   var_marginaly_2d(as.vector(x), var_mean0, var_e, type[1], var_margy0), 
                                                                   var_marginaly_2d(as.vector(x), var_mean1, var_e, type[2], var_margy1)))
  old_initD = initD
  initD = initD[-which(w_initD == 0),]
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
  f0_cand = apply(candidates, 1, FUN = f0)
  f1_cand = apply(candidates, 1, FUN = f1)
  xinitind = which.max(abs(f0_cand - f1_cand)) # same index as max wasserstein distance, but maybe at some point implement that one instead
  x_maxW = candidates[xinitind, ]
  is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), x_maxW))
  if(is_x_max_in_initD){
    N = N - 1
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, initD[-25,], k, mean_beta0, mean_beta1, 
                                                                 var_mean0, var_mean1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p))
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
                                                                 var_mean0, var_mean1, var_e, f0, f1, 
                                                                 type, var_margy0, var_margy1, p))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
    print(i)
  }
  return(list("initD" = old_initD, "addD" = D, "updatedD" = rbind(old_initD, D), "q_initD" = initD))
}

MED_ms_fast_2d_add = function(initD, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                              f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                              N = 11, numCandidates = NULL, K = 10, p = 2, xmin = 0, xmax = 1, seed = 1){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  w_initD = apply(initD, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                   var_marginaly_2d(as.vector(x), var_mean0, var_e, type[1], var_margy0), 
                                                                   var_marginaly_2d(as.vector(x), var_mean1, var_e, type[2], var_margy1)))
  old_initD = initD
  initD = initD[-which(w_initD == 0),]
  # check if any points are equal to point with max wasserstein
  
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
  x_maxW = C[,,1][xinitind, ]
  is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), x_maxW))
  
  if(is_x_max_in_initD){
    # just use x_maxW, but take it out at the end, since it's already included. Note that this means N = N - 1
    Dkplus[1, ] = x_maxW
    initD = initD[-which(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), x_maxW)),]
  } else{
    Dkplus[1, ] = x_maxW # x1, first element of set of design points, D
  }
  
  for(k in 1:(K - 1)){
    
    if(is_x_max_in_initD){
      # just use x_maxW, but take it out at the end, since it's already included. Note that this means N = N - 1
      Dkplus[1, ] = x_maxW
    } else{
      Dkplus[1, ] = x_maxW # x1, first element of set of design points, D
    }
    
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
      
      f_min_candidates = apply(C[ , , j], 1, function(x) f_min_fast_2d(x, rbind(initD, Dkplus[1:(j - 1), , drop = FALSE]), gammas[k], 
                                                                       mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                       f0, f1, type, var_margy0, var_margy1, p))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      Dkplus[j, ] = C[ , , j][chosen_cand, ]
      print(paste("k:", k, " out of ", (K - 1)," --- j:", j, " out of ", N))
    }
    Dk = Dkplus
  }
  return(list("initD" = old_initD, "addD" = Dkplus, "updatedD" = rbind(old_initD, Dkplus), "q_initD" = initD, "candidates" = C))
}











##################
### evaluate D ###
##################

### --- Compute Criteria --- ###

totalPE_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  # if(log_space == FALSE){
  #   numPairs = N * (N - 1) / 2
  #   pairwise_PEs = rep(NA, numPairs)
  #   counter = 1
  #   qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  #   for(i in 1:(N - 1)){
  #     for(j in (i + 1):N){
  #       pairwise_PEs[counter] = qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2))
  #       counter = counter + 1
  #     }
  #   }
  #   return(sum(pairwise_PEs))
  # } else{
    numPairs = N * (N - 1) / 2
    pairwise_terms = rep(NA, numPairs)
    counter = 1
    logqD = apply(D, 1, FUN = function(x) log(q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                                   f0, f1, type, var_margy0, var_margy1, p)))
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_terms[counter] = logqD[i] + logqD[j] - log(sqrt(sum((D[i, ] - D[j, ])^2)))
        counter = counter + 1
      }
    }
    return(exp(logSumExp(pairwise_terms)))
  # }
}


crit_1atatime_2d = function(D, N, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  # if(log_space == FALSE) {
    # numPairs = N * (N - 1) / 2
    # pairwise_PEs = rep(NA, numPairs)
    # counter = 1
    # qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
    # for(i in 1:(N - 1)){
    #   for(j in (i + 1):N){
    #     pairwise_PEs[counter] = (qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2)))^k
    #     counter = counter + 1
    #   }
    # }
    # return((sum(pairwise_PEs))^(1/k))
  # } else{
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
  # }
}

crit_fast_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != dim(D)[1]) stop("N is not the same as length of D")
  # if(log_space == FALSE){
  #   numPairs = N * (N - 1) / 2
  #   pairwise_PEs = rep(NA, numPairs)
  #   counter = 1
  #   qD = apply(D, 1,  FUN = function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p))
  #   for(i in 1:(N - 1)){
  #     for(j in (i + 1):N){
  #       pairwise_PEs[counter] = qD[i] * qD[j] / sqrt(sum((D[i, ] - D[j, ])^2))
  #       counter = counter + 1
  #     }
  #   }
  #   return(max(pairwise_PEs))
  # } else{
    numPairs = N * (N - 1) / 2
    pairwise_terms = rep(NA, numPairs)
    counter = 1
    logqD = apply(D, 1, FUN = function(x) log(q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                                   f0, f1, type, var_margy0, var_margy1, p)))
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_terms[counter] = logqD[i] + logqD[j] - log(sqrt(sum((D[i, ] - D[j, ])^2)))
        counter = counter + 1
      }
    }
    return(max(exp(pairwise_terms)))
  # }
}

model_evidence = function(Y, D, N, mean_beta, var_mean, var_e, type){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  N = length(Y)
  marginaly_mean = X %*% mean_beta
  marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

simulateY = function(D, N, mean_beta, var_mean, var_e, numSims, type = NULL){
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    beta = t(rmvnorm(n = 1, mean = mean_beta, sigma = var_mean))
    for(i in 1:N){
      Y[i, j] = rnorm(n = 1, mean = X[i, ] %*% beta, sd = sqrt(var_e))
    }
  }
  return(Y)
}

calcExpPostProbH_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                               numSims = 100, type = NULL, log_space = TRUE){
  set.seed(123)
  if(log_space == FALSE){
    # --- Y simulated from H0 --- #
    simY = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
      simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      simBF01[j] = simPostH0[j] / simPostH1[j]
    }
    expected_postH0_YH0 = mean(simPostH0)
    expected_postH1_YH0 = mean(simPostH1)
    expected_BF01_YH0 = mean(simBF01)

    # --- Y simulated from H1 --- #
    simY = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
      simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      simBF01[j] = simPostH0[j] / simPostH1[j]
    }
    expected_postH0_YH1 = mean(simPostH0)
    expected_postH1_YH1 = mean(simPostH1)
    expected_BF01_YH1 = mean(simBF01)
  } else{
    # --- Y simulated from H0 --- #
    simY = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
    logsimPostH0 = rep(NA, numSims)
    logsimPostH1 = rep(NA, numSims)
    logsimBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      logsimPostH0[j] = log(simEvidenceH0) - log(simEvidenceH0 + simEvidenceH1)
      logsimPostH1[j] = log(simEvidenceH1) - log(simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      logsimBF01[j] = logsimPostH0[j] - logsimPostH1[j]
    }
    expected_postH0_YH0 = (1 / numSims) * exp(logSumExp(logsimPostH0))
    expected_postH1_YH0 = (1 / numSims) * exp(logSumExp(logsimPostH1))
    expected_BF01_YH0 = (1 / numSims) * exp(logSumExp(logsimBF01))
    
    # --- Y simulated from H1 --- #
    simY = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
    logsimPostH0 = rep(NA, numSims)
    logsimPostH1 = rep(NA, numSims)
    logsimBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      logsimPostH0[j] = log(simEvidenceH0) - log(simEvidenceH0 + simEvidenceH1)
      logsimPostH1[j] = log(simEvidenceH1) - log(simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      logsimBF01[j] = logsimPostH0[j] - logsimPostH1[j]
    }
    expected_postH0_YH1 = (1 / numSims) * exp(logSumExp(logsimPostH0))
    expected_postH1_YH1 = (1 / numSims) * exp(logSumExp(logsimPostH1))
    expected_BF01_YH1 = (1 / numSims) * exp(logSumExp(logsimBF01))
  }
  return(c("expected_postH0_YH0" = expected_postH0_YH0, "expected_postH1_YH0" = expected_postH1_YH0,
           "expected_BF01_YH0" = expected_BF01_YH0, "expected_postH0_YH1" = expected_postH0_YH1,
           "expected_postH1_YH1" = expected_postH1_YH1, "expected_BF01_YH1" = expected_BF01_YH1))
}

# Scaled Prediction Variance (SPV) : N V[y-hat(x_0)] / sigma^2 = N x_0' (X'X)^(-1) x_0

getSPV_2d = function(D, N, type){
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  X = NULL
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
    gridpts = cbind(rep(1, numGridPts), gridpts)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
    gridpts = cbind(rep(1, numGridPts), gridpts[,1], gridpts[,1]^2, gridpts[,2], gridpts[,2]^2)
  }
  XtX = crossprod(X)
  spv = N * apply(gridpts, 1, function(x) x %*% solve(XtX, x))
  spv = sort(spv)
  return(spv)
}


# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise_sqexp2d = function(Xi,Xj, l) exp(-0.5 * sum(((Xi - Xj) / l)^2)) 
C_fn_sqexp2d = function(X, Y, l){
  lenX = dim(X)[1]
  lenY = dim(Y)[1]
  mat = matrix(rep(NA, lenX*lenY), lenX, lenY)
  for(i in 1:lenX){
    for(j in 1:lenY){
      mat[i, j] = C_fn_elementwise_sqexp2d(X[i,], Y[j,], l)
    }
  }
  return(mat)
}

#  here is product exponential
C_fn_elementwise_prodexp2d = function(Xi,Xj, theta) exp(-sum(theta * (Xi - Xj)^2)) 
C_fn_prodexp2d = function(X, Y, theta){
  lenX = dim(X)[1]
  lenY = dim(Y)[1]
  mat = matrix(rep(NA, lenX*lenY), lenX, lenY)
  for(i in 1:lenX){
    for(j in 1:lenY){
      mat[i, j] = C_fn_elementwise_prodexp(X[i,], Y[j,], l)
    }
  }
  return(mat)
}

# doesn't work... det(solve(K_X)) = 0
getSPV_GASP_2d = function(D, N, hyperparam = NULL, type, covtype = 2){
  D = unique(D)
  if(covtype == 1){ # doesn't work yet
    C_fn_elementwise = C_fn_elementwise_prodexp2d
    C_fn = C_fn_prodexp2d
  } else if(covtype == 2){
    C_fn_elementwise = C_fn_elementwise_sqexp2d
    C_fn = C_fn_sqexp2d
  }
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  X = NULL
  if(type[2] == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
    gridpts = cbind(rep(1, numGridPts), gridpts)
    if(is.null(hyperparam)) hyperparam = rep(0.1, 3)
  }
  if(type[2] == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
    gridpts = cbind(rep(1, numGridPts), gridpts[,1], gridpts[,1]^2, gridpts[,2], gridpts[,2]^2)
    if(is.null(hyperparam)) hyperparam = rep(0.1, 5)
  }
  K_X = C_fn(X, X, hyperparam)
  ones = rep(1, N)
  spv = rep(NA, numGridPts)
  for(i in 1:numGridPts){
    newx = gridpts[i,]
    k_newx = apply(X, 1, function(x) C_fn_elementwise(x, newx, hyperparam))
    term2 = t(k_newx) %*% solve(K_X, k_newx)
    term3 = (1 - t(ones) %*% solve(K_X, k_newx))^2 / (t(ones) %*% solve(K_X, ones))
    spv[i] = 1 - term2 + term3
  }
  return(spv)
}




# RMSE for surface h1
RMSE = function(pred, obs){
  sqrt(mean((pred - obs)^2))
}

getRMSE_2d = function(D, N, f_ell){
  # create hypothetical response variable based on test function, here we use f_ell
  fake_y = apply(D, 1, function(x) f_ell(x))
  # analyze design using 2nd order polynomial
  fake_data = as.data.frame(cbind(D, fake_y))
  colnames(fake_data) = c("x1", "x2", "y")
  polymodel_deg1 = lm(fake_y ~ 1 + x1 + x2, data = fake_data)
  polymodel_deg2 = lm(fake_y ~ 1 + x1 + x2 + I(x1^2) + I(x2^2), data = fake_data)
  polymodel_deg3 = lm(fake_y ~ 1 + x1 + x2 + I(x1^2) + I(x2^2) + I(x1^3) + I(x2^3), data = fake_data)
  # predict responses across grid
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  new_data = as.data.frame(gridpts)
  colnames(new_data) = c("x1", "x2")
  polymodel_deg1_predy = predict(polymodel_deg1, newdata = new_data)
  polymodel_deg2_predy = predict(polymodel_deg2, newdata = new_data)
  polymodel_deg3_predy = predict(polymodel_deg3, newdata = new_data)
  # calculate RMSE
  obs = apply(gridpts, 1, function(x) f_ell(x))
  polymodel_deg1_RMSE = RMSE(polymodel_deg1_predy, obs)
  polymodel_deg2_RMSE = RMSE(polymodel_deg2_predy, obs)
  polymodel_deg3_RMSE = RMSE(polymodel_deg3_predy, obs)
  return(list("RMSE_deg1" = polymodel_deg1_RMSE, "RMSE_deg2" = polymodel_deg2_RMSE, "RMSE_deg3" = polymodel_deg3_RMSE))
}

