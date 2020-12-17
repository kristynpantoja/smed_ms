

#############
#############
### D > 1 ###
#############
#############

###############################################################
### with data (uses posterior predictive distribution of y) ###
###############################################################

f_min_vs = function(candidate, D, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, k, alpha, buffer){
  result = q_vs(candidate, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, 
                alpha, buffer)^k * 
    sum(apply(D, 1, function(x_i) (q_vs(x_i, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, 
                                        alpha, buffer) / sqrt(sum((x_i - candidate)^2)))^k))
  return(result)
}

# multidimensional case
# add_MED_ms_oneatatime_data_multidim
add_MMEDvs_oneatatime = function(initD, inity, mean_beta_full, beta_true = NULL, indices_true = NULL, 
                                 indices0, indices1, mean_beta0 = NULL, mean_beta1 = NULL, 
                                 var_e = 1, var_beta = NULL, var_beta0 = NULL, var_beta1 = NULL,
                                 N2 = 11, xmin = -1, xmax = 1, numCandidates = 10^5, k = 4, p = 1, 
                                 alpha = 1, buffer = 0, candidates = NULL,
                                 wasserstein0 = 1, genCandidates = 1, log_space = FALSE){
  initN = dim(initD)[1]
  if(length(inity) != initN) stop("length of inity does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta
  postvar0 = postvar(initD[ , indices0], initN, var_e, var_beta0, type = NULL)
  postmean0 = postmean(inity, initD[ , indices0], initN, mean_beta0, var_beta0, var_e, type = NULL)
  postvar1 = postvar(initD[ , indices1], initN, var_e, var_beta1, type = NULL)
  postmean1 = postmean(inity, initD[ , indices1], initN, mean_beta1, var_beta1, var_e, type = NULL)
  
  # do something about bad initial data points that will make the TPE explode
  if(wasserstein0 == 1){
    w_initD = apply(initD, 1, FUN = function(x) Wasserstein_vs(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e))
    if(length(which(w_initD == 0)) != 0){
      initD = initD[-which(w_initD == 0), ]
      y = inity[-which(w_initD == 0)]
      w_initD = w_initD[-which(w_initD == 0)]
      postvar0 = postvar(initD[ , indices0], initN, var_e, var_beta0, type = NULL)
      postmean0 = postmean(inity, initD[ , indices0], initN, mean_beta0, var_beta0, var_e, type = NULL)
      postvar1 = postvar(initD[ , indices1], initN, var_e, var_beta1, type = NULL)
      postmean1 = postmean(inity, initD[ , indices1], initN, mean_beta1, var_beta1, var_e, type = NULL)
    }
  }
  if(wasserstein0 == 2){
    if(buffer == 0) warning("Buffer = 0, but wasserstein0 = 2; will assign Buffer = 1e-4")
    buffer == 1e-4
  }
  
  # other variables and checks
  ttlN = initN + N2
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    # if(genCandidates == 1) candidates = expand.grid(seq(from = xmin, to = xmax, length.out = floor(numCandidates^(1 / length(mean_beta_full)))))
    candidates = matrix(runif(n = pfull * numCandidates, min = xmin, max = xmax), 
                        nrow = numCandidates, ncol = pfull)
  }
  
  # -- Initialize 1st additional design point-- #
  D = matrix(rep(NA, N2 * pfull), N2, pfull)
  w_cand = apply(candidates, 1, FUN = function(x) Wasserstein_vs(x, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e))
  xinitind = which.max(w_cand)
  xmaxW = candidates[xinitind, ]
  is_x_max_in_initD = any(apply(initD, 1, function(x, want) isTRUE(all.equal(x, want)), xmaxW)) # give tolerance?
  if(is_x_max_in_initD){
    # stop(" where - 1 happens? why 1 off?")
    # ????? double check at some point
    # N2 = N2 - 1
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_vs(x, initD, indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, k, alpha, buffer)) # why was this running before changing sapply to apply?
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[1, ] = xnew
  } else{
    D[1, ] = xmaxW
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_vs(x, rbind(initD, D[1:(i - 1), , drop = FALSE]), indices0, indices1, postmean0, postmean1, postvar0, postvar1, var_e, p, k, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
    # print(i)
  }
  
  return(list("initD" = initD, "addD" = D, "D" = rbind(initD, D)))
}


add_MMEDvs = function(initD, inity, mean_beta_full, beta_true = NULL, indices_true = NULL, 
                      indices0, indices1, mean_beta0 = NULL, mean_beta1 = NULL, 
                      var_e = 1, var_beta = NULL, var_beta0 = NULL, var_beta1 = NULL,
                      N2 = 11, xmin = 0, xmax = 1, numCandidates = 10^5, k = 4, p = 1, 
                      alpha = 1, buffer = 0, candidates = NULL,
                      wasserstein0 = 1, genCandidates = 1, log_space = FALSE, 
                      algorithm = 1){
  if(algorithm == 1){
    add_MMEDvs_oneatatime(initD, inity, mean_beta_full, beta_true, indices_true, 
                          indices0, indices1, mean_beta0, mean_beta1, 
                          var_e, var_beta, var_beta0, var_beta1,
                          N2, xmin, xmax, numCandidates, k, p, 
                          alpha, buffer, candidates,
                          wasserstein0, genCandidates, log_space)
  }
}