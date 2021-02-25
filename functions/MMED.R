# when no input data is given, use MED on the marginal distribution of y 
#   (prior predictive)
objective_mmed = function(
  candidate, D, k, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, 
  f0, f1, type, var_margy0, var_margy1, p, alpha, log_space = FALSE
){
  if(log_space == FALSE) {
    potential.energy = function(x){
      (q_mmed(x, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, 
              f0, f1, type, var_margy0, var_margy1, p, alpha) / 
         sqrt((x - candidate)^2))^k
    }
    
    result = q_mmed(
      candidate, beta.mean0, beta.mean1, beta.var0, beta.var1, 
      error.var, f0, f1, type, var_margy0, var_margy1, p, alpha)^k * 
      sum(sapply(D, FUN = potential.energy))
    return(result)
  } else{
    # if has logSumExp library
    log.potential.energy = function(x){
      k * 
        log(q(x, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, 
              f1, type, var_margy0, var_margy1, p, alpha)) - 
        (k / 2) * log((x - candidate)^2)
    }
    candidate.logpe = k * 
      log(q(candidate, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, 
            f0, f1, type, var_margy0, var_margy1, p, alpha))
    D.logpe = sapply(D, FUN = log.potential.energy)
    result = exp(logSumExp(candidate.logpe + D.logpe))
    return(result)
  }
}

MMED = function(
  beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, 
  f0 = NULL, f1 = NULL, type = NULL, N = 11, numCandidates = 10^5, k = 4, 
  xmin = 0, xmax = 1, p = 2, alpha = NULL, genCandidates = 1, initialpt = 1, 
  var_margy0 = NULL, var_margy1 = NULL, 
  log_space = FALSE
){
  # some error checking, first:
  if(is.null(type) & 
     is.null(f0) & 
     is.null(f1) & 
     is.null(var_margy0) 
     & is.null(var_margy0)) stop("must specify model type and/or model")
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
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  old_candidates = candidates
  # -- Initialize 1st Design Point in D -- #
  D = rep(NA, N)
  if(initialpt == 1){
    wassfn = function(x){
      mu1.temp = f0(x) # mean of marginal dist of y | H0
      mu2.temp = f1(x) # mean of marginal dist of y | H1
      var1.temp = var_marginaly(x, V0, sigmasq, type = type01[1]) # variance of marginal dist of y | H0
      var2.temp = var_marginaly(x, V1, sigmasq, type = type01[2]) # variance of marginal dist of y | H1
      WN(mu1.temp, mu2.temp, var1.temp, var2.temp)
    }
    w_vec = sapply(candidates, wassfn)
    xopt.idx = which.max(w_vec)
    D[1] = candidates[xopt.idx]
  }
  
  if(N > 1){
    for(i in 2:N){
      # Find f_opt: minimum of f_min
      # calculate cumulative TPE
      f_min_candidates = sapply(
        candidates, 
        function(x) objective_mmed(
          x, D[1:(i - 1)], k, beta.mean0, beta.mean1, beta.var0, beta.var1, 
          error.var, f0, f1, type, var_margy0, var_margy1, p, alpha, 
          log_space))
      f_opt = which.min(f_min_candidates)
      D[i] = candidates[f_opt]
    }
  }
  return(D)
}
