# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("variance_marginal_y.R")

##########
### 2D ###
##########

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

