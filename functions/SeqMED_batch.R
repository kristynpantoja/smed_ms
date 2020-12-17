# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("variance_marginal_y.R")
# require("construct_design_matrix.R")
# require("posterior_variance.R")
# require("posterior_mean.R")

##########
##########
### 1D ###
##########
##########



######################################################
### without data (uses marginal distribution of y) ###
######################################################
# MMED_batch = function(
#   initD, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e,
#   f0 = NULL, f1 = NULL, type = NULL, N2 = 11, numCandidates = 10^5, k = 4,
#   xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0,
#   wasserstein0 = 1, genCandidates = 1, candidates = NULL,
#   var_margy0 = NULL, var_margy1 = NULL, log_space = FALSE
# ){
#   # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
# 
#   # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
#   old_initD = initD
#   if(wasserstein0 == 1){
#     w_initD = sapply(initD, FUN = function(x) WN(f0(x), f1(x),
#                                                  var_marginaly(x, var_beta0, var_e, type[1], var_margy0),
#                                                  var_marginaly(x, var_beta1, var_e, type[2], var_margy1)))
#     if(length(which(w_initD == 0)) != 0){
#       initD = initD[-which(w_initD == 0)]
#       w_initD = w_initD[-which(w_initD == 0)]
#     }
#   }
#   if(wasserstein0 == 2){
#     if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 0.0001")
#     buffer == 0.0001
#   }
# 
#   # other variables and checks
#   initN = length(initD)
#   ttlN = initN + N2
#   if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
# 
#   # Create hypothesized models
#   if(is.null(f0)){
#     if(type[1] == 1) f0 = function(x) mean_beta0 * x
#     else if(type[1] == 2) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x
#     else if(type[1] == 3) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x + mean_beta0[3] * x^2
#     else stop("type[1] is invalid and f0 is not provided")
#   }
#   if(is.null(f1)){
#     if(type[2] == 1) f1 = function(x) mean_beta1 * x
#     else if(type[2] == 2) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x
#     else if(type[2] == 3) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2
#     else stop("type[2] is invalid and f1 is not provided")
#   }
# 
#   # -- Generate Candidate Points -- #
#   if(is.null(candidates)){
#     if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
#     if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
#   }
# 
#   # -- Initialize 1st additional design point-- #
#   D = rep(NA, N2)
#   optimal_q = optimize(function(x) q_mmed(x, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1,
#                                      type, var_margy0, var_margy1, p, alpha, buffer), interval = c(xmin, xmax))$minimum
#   xopt = optimal_q
#   is_x_max_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
#   if(is_x_max_in_initD){
#     # Find f_opt: minimum of f_min
#     f_min_candidates = sapply(candidates, function(x) f_min_mmed(x, initD, k, mean_beta0, mean_beta1,
#                                                                  var_beta0, var_beta1, var_e, f0, f1,
#                                                                  type, var_margy0, var_margy1, p, alpha, buffer))
#     f_opt = which.min(f_min_candidates)
#     xnew = candidates[f_opt]
#     # Update set of design points (D) and plot new point
#     D[1] = xnew
#   } else{
#     D[1] = xopt
#   }
# 
#   for(i in 2:N2){
#     # Find f_opt: minimum of f_min
#     f_min_candidates = sapply(candidates, function(x) f_min_mmed(x, c(initD, D[1:(i - 1)]), k, mean_beta0, mean_beta1,
#                                                                    var_beta0, var_beta1, var_e, f0, f1,
#                                                                    type, var_margy0, var_margy1, p, alpha, buffer))
#     f_opt = which.min(f_min_candidates)
#     xnew = candidates[f_opt]
#     # Update set of design points (D) and plot new point
#     D[i] = xnew
#   }
#   return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD))
# }



###############################################################
### with data (uses posterior predictive distribution of y) ###
###############################################################

f_min_seqmed = function(candidate, D, postmean0, postmean1, postvar0, postvar1, var_e, type, p, k, alpha, buffer){
  result = q_seqmed(candidate, postmean0, postmean1, postvar0, postvar1, var_e, type, p, 
                  alpha, buffer)^k * 
    sum(sapply(D, function(x_i) (q_seqmed(x_i, postmean0, postmean1, postvar0, postvar1, var_e, type, p, 
                                        alpha, buffer) / sqrt((x_i - candidate)^2))^k))
  return(result)
}

# add_MED_ms_oneatatime_data
SeqMED_batch= function(
  initD, y, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
  f0 = NULL, f1 = NULL, type = NULL, N2 = 11, numCandidates = 10^5, k = 4, 
  xmin = -1, xmax = 1, p = 1, alpha = NULL, buffer = 0, 
  wasserstein0 = 1, genCandidates = 1, candidates = NULL, 
  var_margy0 = NULL, var_margy1 = NULL, log_space = FALSE
){
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
    w_initD = sapply(initD, FUN = function(x) WNlm(x, postmean0, postmean1, postvar0, postvar1, var_e, type))
    if(length(which(w_initD == 0)) != 0){
      initD = initD[-which(w_initD == 0)]
      y = y[-which(w_initD == 0)]
      w_initD = w_initD[-which(w_initD == 0)]
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
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & 
     is.null(var_margy0)) stop("must specify model type and/or model")
  
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
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  optimal_q = optimize(function(x) q_seqmed(x, postmean0, postmean1, postvar0, postvar1, var_e, type, p,
                                           alpha, buffer), interval = c(xmin, xmax))$minimum
  xopt = optimal_q
  is_x_max_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_seqmed(x, initD, postmean0, postmean1, postvar0, postvar1, var_e, type, p, k, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
  } else{
    D[1] = xopt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = sapply(candidates, function(x) f_min_seqmed(x, c(initD, D[1:(i - 1)]), postmean0, postmean1, postvar0, postvar1, var_e, type, p, k, alpha, buffer))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt]
      # Update set of design points (D) and plot new point
      D[i] = xnew
    }
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD))
}



















# 
# ##########
# ### 2D ###
# ##########
# 
# # hasn't been tested
# add_MED_ms_oneatatime_2d = function(initD, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
#                                         f0 = NULL, f1 = NULL, type = NULL, N = 11, numCandidates = 10^5, k = 4, 
#                                         xmin = 0, xmax = 1, p = 3, alpha = NULL, buffer = 0, y = NULL, 
#                                         log_space = FALSE, genCandidates = 1, wasserstein0 = 1, 
#                                         var_margy0 = NULL, var_margy1 = NULL){
#   # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
#   # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
#   old_initD = initD
#   
#   if(wasserstein0 == 1){
#     w_initD = apply(initD, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
#                                                                      var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
#                                                                      var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
#     initD = initD[-which(w_initD == 0), ]
#     w_initD = w_initD[-which(w_initD == 0), ]
#   }
#   if(wasserstein0 == 2){
#     if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 0.0001")
#     buffer == 0.0001
#   }
#   # other variables and checks
#   initN = dim(initD)[1]
#   ttlN = initN + N
#   if(dim(initD)[2] != 2) stop("initD is not a matrix of size initN x 2")
#   if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
#   
#   # Create hypothesized models
#   if(is.null(f0)){
#     if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
#     else if(type[1] == 5) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[1]^2 + mean_beta0[4] * x[2] + mean_beta0[3] * x[2]^2
#     else stop("type[1] is invalid and f0 is not provided")
#   }
#   if(is.null(f1)){
#     if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
#     else if(type[1] == 5) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[1]^2 + mean_beta1[4] * x[2] + mean_beta1[3] * x[2]^2
#     else stop("type[2] is invalid and f1 is not provided")
#   }
#   
#   # -- Generate Candidate Points -- #
#   
#   if(length(numCandidates) == 1) numCandidates = c(floor(sqrt(numCandidates)), floor(sqrt(numCandidates)))
#   # sometimes one factor is more important, so gets more candidates, than the other factor.
#   if(genCandidates == 1){
#     candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
#     candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
#     candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
#                        rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
#   }
#   if(genCandidates == 2){
#     candidates_x1 = sort(runif(numCandidates[1], min = xmin, max = xmax))
#     candidates_x2 = sort(runif(numCandidates[2], min = xmin, max = xmax))
#     candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
#                        rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
#   }
#   
#   # -- Initialize 1st Design Point in D -- #
#   D = matrix(rep(NA, N * 2), N, 2)
#   w_cand = apply(candidates, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
#                                                                        var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
#                                                                        var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
#   xinitind = which.max(w_cand)
#   xmaxW = candidates[xinitind, ]
#   is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), xmaxW)) # give tolerance?
#   if(is_x_max_in_initD){
#     N = N - 1
#     # Find f_opt: minimum of f_min
#     f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, initD, k, mean_beta0, mean_beta1, 
#                                                                  var_beta0, var_beta1, var_e, f0, f1, 
#                                                                  type, var_margy0, var_margy1, p, alpha, buffer))
#     f_opt = which.min(f_min_candidates)
#     xnew = candidates[f_opt, ]
#     # Update set of design points (D) and plot new point
#     D[1, ] = xnew
#   } else{
#     D[1, ] = x_maxW
#   }
#   
#   for(i in 2:N){
#     # Find f_opt: minimum of f_min
#     f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, rbind(initD, D[1:(i - 1), , drop = FALSE]), k, mean_beta0, mean_beta1, 
#                                                                  var_beta0, var_beta1, var_e, f0, f1, 
#                                                                  type, var_margy0, var_margy1, p, alpha, buffer))
#     f_opt = which.min(f_min_candidates)
#     xnew = candidates[f_opt, ]
#     # Update set of design points (D) and plot new point
#     D[i, ] = xnew
#     print(i)
#   }
#   return(list("initD" = old_initD, "addD" = D, "updatedD" = rbind(old_initD, D), "q_initD" = initD))
# }
# 
# 
# 
# # require("wasserstein_distance.R")
# # require("charge_function_q.R")
# # require("variance_marginal_y.R")
# 
# ##########
# ### 2D ###
# ##########
# 
# # not tested yet
# augment_MED_ms_fast_2d = function(initD, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
#                                   f0 = NULL, f1 = NULL, type = NULL, N = 11, K = 10, 
#                                   xmin = 0, xmax = 1, p = 2, alpha = NULL, buffer = 0, numCandidates = NULL,
#                                   var_margy0 = NULL, var_margy1 = NULL, seed = 1){
#   # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
#   
#   # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
#   w_initD = apply(initD, 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
#                                                                    var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
#                                                                    var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
#   old_initD = initD
#   initD = initD[-which(w_initD == 0),]
#   # check if any points are equal to point with max wasserstein
#   
#   # other variables and checks
#   initN = dim(initD)[1]
#   ttlN = initN + N
#   if(dim(initD)[2] != 2) stop("initD is not a matrix of size initN x 2")
#   if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
#   
#   # Create hypothesized models
#   if(is.null(f0)){
#     if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
#     else if(type[1] == 5) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[1]^2 + mean_beta0[4] * x[2] + mean_beta0[3] * x[2]^2
#     else stop("type[1] is invalid and f0 is not provided")
#   }
#   if(is.null(f1)){
#     if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
#     else if(type[1] == 5) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[1]^2 + mean_beta1[4] * x[2] + mean_beta1[3] * x[2]^2
#     else stop("type[2] is invalid and f1 is not provided")
#   }
#   
#   # -- Make D_1, space-filling design -- #
#   # check that n >= 3
#   if(N < 3) stop("not enough samples - need at least 3.")
#   
#   # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
#   sqrtNumCand = NULL
#   if(is.null(numCandidates)){
#     # make the number of candidates approximately equal to N (must be an integer square root, to make a grid)
#     sqrtNumCand = floor(sqrt(N))
#   }
#   if(length(numCandidates) == 1){ # 
#     sqrtNumCand = floor(sqrt(numCandidates)) # to be able to make a grid of candidates
#   }
#   numCandidates = c(sqrtNumCand, sqrtNumCand)
#   numCandidatesTtl = sqrtNumCand^2
#   candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
#   candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
#   C1 = cbind(rep(candidates_x1, each = numCandidates[2]), 
#              rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
#   spots_left = N - dim(C1)[1]
#   if(spots_left != 0){
#     set.seed(seed)
#     candidates_x1_leftover = runif(spots_left, xmin, xmax)
#     candidates_x2_leftover = runif(spots_left, xmin, xmax)
#     D1 =  rbind(C1, cbind(candidates_x1_leftover, candidates_x2_leftover))
#   } else{
#     D1 = C1
#   }
#   
#   ## calculate Wassersteins for each design point
#   #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
#   # if do this instead of calculating q for each pair each time, may save time and space
#   
#   # -- If K = 1, return the space-filling design -- #
#   if(K == 1){
#     D = D1
#     C = C1
#     return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
#   }
#   
#   # -- If K > 1, choose new design -- #
#   Dkplus = matrix(rep(NA, N), N, 2)
#   Dk = D1
#   #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
#   gammas = c(1:(K - 1)) / (K - 1)
#   # -- the final step should be gamma = 1 because then we optimize the correct criterion
#   # save candidates for each K
#   C <- array(NA, dim = c(numCandidatesTtl * K, 2, N))
#   for (j in 1:N){
#     C[ , , j] = rbind(C1, matrix(rep(NA, 2 * numCandidatesTtl * (K - 1)), numCandidatesTtl * (K - 1), 2))
#   }
#   
#   # at index k, determine the next design k + 1
#   
#   ## For j = 1, i.e. 1st point in design k + 1: it's always the same point for all k ##
#   # criterion to choose first candidate from candidate set: 
#   # the point at which f1 and f2 are most different
#   w_evals = apply(C[,,1], 1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
#                                                                     var_marginaly_2d(as.vector(x), var_beta0, var_e, type[1], var_margy0), 
#                                                                     var_marginaly_2d(as.vector(x), var_beta1, var_e, type[2], var_margy1)))
#   # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
#   xinitind = which.max(w_evals)
#   x_maxW = C[,,1][xinitind, ]
#   is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), x_maxW))
#   
#   if(is_x_max_in_initD){
#     # just use x_maxW, but take it out at the end, since it's already included. Note that this means N = N - 1
#     Dkplus[1, ] = x_maxW
#     initD = initD[-which(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), x_maxW)),]
#   } else{
#     Dkplus[1, ] = x_maxW # x1, first element of set of design points, D
#   }
#   
#   for(k in 1:(K - 1)){
#     
#     if(is_x_max_in_initD){
#       # just use x_maxW, but take it out at the end, since it's already included. Note that this means N = N - 1
#       Dkplus[1, ] = x_maxW
#     } else{
#       Dkplus[1, ] = x_maxW # x1, first element of set of design points, D
#     }
#     
#     # for j = 2:N
#     for(j in 2:N){
#       # get candidates in neighborhood L_jk = (lower, upper)
#       diff = (Dk[-j, ] - matrix(rep(Dk[j, ], N - 1), N - 1, 2, byrow = TRUE))^2
#       Rjk = min(sqrt(rowSums(diff)))
#       lower_x1 = max(Dk[j, 1] - Rjk, 0); lower_x2 = max(Dk[j, 2] - Rjk, 0)
#       upper_x1 = min(Dk[j, 1] + Rjk, 1); upper_x2 = min(Dk[j, 2] + Rjk, 1)
#       tildeDj_kplus_x1 = seq(from = lower_x1, to = upper_x1, length.out = numCandidates[1])
#       tildeDj_kplus_x2 = seq(from = lower_x2, to = upper_x2, length.out = numCandidates[2])
#       tildeDj_kplus = cbind(rep(tildeDj_kplus_x1, each = numCandidates[2]), 
#                             rep(tildeDj_kplus_x2, times = numCandidates[1]))
#       C[(k * numCandidatesTtl + 1):((k + 1) * numCandidatesTtl) , , j] = tildeDj_kplus # This is now C_j^{k+1}
#       
#       f_min_candidates = apply(C[ , , j], 1, function(x) f_min_fast_2d(x, rbind(initD, Dkplus[1:(j - 1), , drop = FALSE]), gammas[k], 
#                                                                        mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
#                                                                        f0, f1, type, var_margy0, var_margy1, p, alpha, buffer))
#       #choose that which has largest evaluation of criterion
#       chosen_cand = which.min(f_min_candidates)
#       Dkplus[j, ] = C[ , , j][chosen_cand, ]
#       print(paste("k:", k, " out of ", (K - 1)," --- j:", j, " out of ", N))
#     }
#     Dk = Dkplus
#   }
#   return(list("initD" = old_initD, "addD" = Dkplus, "updatedD" = rbind(old_initD, Dkplus), "q_initD" = initD, "candidates" = C))
# }