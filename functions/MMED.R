# when no input data is given, use MED on the marginal distribution of y 
#   (prior predictive)
objective_mmed = function(
  candidate, D, k, model0, model1, error.var, 
  var_margy0, var_margy1, p, alpha #, log_space = FALSE
){
  # if(log_space == FALSE) {
    potential.energy = function(x){
      (q_mmed(x, model0, model1, error.var, var_margy0, var_margy1, p, alpha) / 
         sqrt((x - candidate)^2))^k
    }
    
    result = q_mmed(
      candidate, model0, model1, 
      error.var, var_margy0, var_margy1, p, alpha)^k * 
      sum(sapply(D, FUN = potential.energy))
    return(result)
}

MMED = function(
  model0, model1, error.var, N = 11, numCandidates = 10^5, k = 4, 
  xmin = 0, xmax = 1, p = 1, alpha = NULL, genCandidates = 1, initialpt = 1, 
  var_margy0 = NULL, var_margy1 = NULL #, log_space = FALSE
){
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  old_candidates = candidates
  
  # check if var_margy0 and var_margy1 are given
  if(is.null(var_margy0) | is.null(var_margy1)){
    # assuming intercept
    var_marg = function(x, model){
      X = model$designMat(x)
      pow = c(1, 2 * c(1:(ncol(X) - 1)))
      return(X^pow %*% diag(model$beta.var))
    }
    if(is.null(var_margy0)) var_margy0 = function(x) var_marg(x, model0)
    if(is.null(var_margy1)) var_margy1 = function(x) var_marg(x, model1)
  }
  
  # -- Initialize 1st Design Point in D -- #
  D = rep(NA, N)
  if(initialpt == 1){
    # assuming intercept term
    wassfn = function(x){
      mu1.temp = model0$designMat(x) %*% model0$beta.mean # mean of marginal dist of y | H0
      mu2.temp = model1$designMat(x)  %*% model1$beta.mean # mean of marginal dist of y | H1
      var1.temp = var_margy0(x) # variance of marginal dist of y | H0
      var2.temp = var_margy1(x) # variance of marginal dist of y | H1
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
          x, D[1:(i - 1)], k, model0, model1, 
          error.var, var_margy0, var_margy1, p, alpha))
      f_opt = which.min(f_min_candidates)
      D[i] = candidates[f_opt]
    }
  }
  return(D)
}
