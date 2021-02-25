objective_seqmed = function(
  candidate, D, postmean0, postmean1, postvar0, postvar1, error.var, type, 
  p = 1, k = 4, alpha = 1){
  result = q_seqmed(
    candidate, postmean0, postmean1, postvar0, postvar1, error.var, type, p, 
    alpha)^k * 
    sum(sapply(
      D, 
      function(x_i) (q_seqmed(x_i, postmean0, postmean1, postvar0, postvar1, 
                              error.var, type, p, alpha) / 
                       sqrt((x_i - candidate)^2))^k))
  return(result)
}

SeqMED_batch = function(
  initD, y, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, 
  f0 = NULL, f1 = NULL, type = NULL, N2 = 11, numCandidates = 10^5, k = 4, 
  xmin = -1, xmax = 1, p = 1, alpha = 1, genCandidates = 1, candidates = NULL, 
  batch.idx = 1
){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta
  postvar0 = postvar(initD, initN, error.var, beta.var0, type[1])
  postmean0 = postmean(y, initD, initN, beta.mean0, beta.var0, error.var, type[1])
  postvar1 = postvar(initD, initN, error.var, beta.var1, type[2])
  postmean1 = postmean(y, initD, initN, beta.mean1, beta.var1, error.var, type[2])
  
  w_initD = sapply(initD, FUN = function(x) WNlm(
    x, postmean0, postmean1, postvar0, postvar1, error.var, type))
  # take out the values with wasserstein distance equal to 0
  if(length(which(w_initD == 0)) != 0){
    initD = initD[-which(w_initD == 0)]
    y = y[-which(w_initD == 0)]
    w_initD = w_initD[-which(w_initD == 0)]
    postvar0 = postvar(initD, initN, error.var, beta.var0, type[1])
    postmean0 = postmean(y, initD, initN, beta.mean0, beta.var0, error.var, type[1])
    postvar1 = postvar(initD, initN, error.var, beta.var1, type[2])
    postmean1 = postmean(y, initD, initN, beta.mean1, beta.var1, error.var, type[2])
  }
  
  # other variables and checks
  initN = length(initD)
  ttlN = initN + N2
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 1) f0 = function(x) beta.mean0 * x
    else if(type[1] == 2) f0 = function(x) beta.mean0[1] + beta.mean0[2] * x
    else if(type[1] == 3) f0 = function(x) beta.mean0[1] + beta.mean0[2] * x + beta.mean0[3] * x^2
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[2] == 1) f1 = function(x) beta.mean1 * x
    else if(type[2] == 2) f1 = function(x) beta.mean1[1] + beta.mean1[2] * x
    else if(type[2] == 3) f1 = function(x) beta.mean1[1] + beta.mean1[2] * x + beta.mean1[3] * x^2
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  D = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, function(x) WNlm(
      x, postmean0, postmean1, postvar0, postvar1, error.var, type))
    w_opt = which.max(w_candidates)
    xopt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))
  } else{
    is_x_max_in_initD = TRUE
  }
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) objective_seqmed(x, initD, postmean0, postmean1, postvar0, 
                               postvar1, error.var, type, p, k, alpha))
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
      f_min_candidates = sapply(
        candidates, 
        function(x) objective_seqmed(x, c(initD, D[1:(i - 1)]), postmean0, 
                                 postmean1, postvar0, postvar1, error.var, type, 
                                 p, k, alpha))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt]
      # Update set of design points (D) and plot new point
      D[i] = xnew
    }
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD))
}
