# Wasserstein distance betwen two normals
Wasserstein_distance = function(mu1, mu2, var1, var2, dim = 1){
  # Normal(mu1, var1)
  # Normal(mu2, var2)
  # dim:
  #   1 is for univariate case
  #   > 1 is for multivariate case
  if(dim == 1){
    wass = sqrt((mu1 - mu2)^2 + (sqrt(var1) - sqrt(var2))^2)
  } else{
    if(dim > 1){
      sqrt_var2 = sqrtm(var2)
      wass = sqrt(crossprod(mu1 - mu2) + sum(diag(var1 + var2 - 2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2))))
    } else{
      stop("invalid dim")
    }
  }
  return(as.numeric(wass))
}

# Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly = function(x, var_mean, var_e, type, var_margy){
  # type:
  #   1 for linear model without slope
  #   2 for linear model with slope
  #   3 for quadratic model with slope
  if(!is.null(type)){
    if(type == 1) var_e + x^2 * var_mean
    else if(type == 2) var_mean[1] + x^2 * var_mean[2] + var_e
    else if(type == 3) var_mean[1] + x^2 * var_mean[2] + x^4 * var_mean[3] + var_e
    else stop(paste("invalid type given : ", type))
  } else{
    var_margy(x = x, var_mean = var_mean, var_e = var_e)
  }
  
}

# charge function evaluated at x
q = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(1.0 / Wass_dist)
}



###################
### one at time ###
###################

f_min = function(candidate, D, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                 f0, f1, type, var_margy0, var_margy1, p, log_space = FALSE){
  if(log_space == FALSE) {
    result = q(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)^k * sum(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p) / sqrt((x_i - candidate)^2))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) + k * log(q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) - k * log(sqrt((x_i - candidate)^2)))
    result = exp(logSumExp(terms))
    return(result)
  }
}


MED_ms = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                   f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                   N = 11, numCandidates = 10^5, k = 4, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
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
    f_min_candidates = sapply(candidates, function(x) f_min(x, D[1:(i - 1)], k, mean_beta0, mean_beta1, 
                                                            var_mean0, var_mean1, var_e, f0, f1, 
                                                            type, var_margy0, var_margy1, p, log_space))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
  }
  return(D)
}



##################
###    fast    ###
##################

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



##################
### evaluate D ###
##################

### --- Posterior Variance --- ###

postvar = function(D, N, var_e, var_mean, type, diagPrior = TRUE){
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * (1 / var_mean)))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * diag(1/diag(var_mean))))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  }
}

# scaled prediction variance
getSPV = function(X, N){
  XtX = t(X) %*% X
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    SPV = sapply(X, FUN = function(x) N * t(as.matrix(x)) %*% solve(XtX, x))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    XtX = t(X) %*% X
    SPV = apply(X, 1, FUN = function(x) N * t(as.matrix(x)) %*% solve(XtX, x))
  }
  return(SPV)
}

# variance dispersion graphs (VDGs) to evaluate the prediction variance properties for experimental designs.
# allows the practitioner to gain an impression regarding the stability of the scaled prediction variance

# the ratio of the prediction variance to the error variance
# is not a function of the error variance. This ratio, called the relative variance of prediction, depends only on
# the design and the factor setting and can be calculated before acquiring the data.
# The prediction variance profile plots the relative variance of prediction as a function of 
# each factor at fixed values of the other factors 



### --- Compute Criterion --- ###

totalPE = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / sqrt((D[i] - D[j])^2)
      counter = counter + 1
    }
  }
  return(sum(pairwise_PEs))
}


crit_1atatime = function(D, N, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p, log_space = TRUE){
  if(N != length(D)) stop("N is not the same as length of D")
  if(log_space == FALSE) {
    numPairs = N * (N - 1) / 2
    pairwise_PEs = rep(NA, numPairs)
    counter = 1
    qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p), D)
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_PEs[counter] = (qD[i] * qD[j] / sqrt((D[i] - D[j])^2))^k
        counter = counter + 1
      }
    }
    return((sum(pairwise_PEs))^(1/k))
  } else{
    if(N != length(D)) stop("N is not the same as length of D")
    numPairs = N * (N - 1) / 2
    pairwise_terms = rep(NA, numPairs)
    counter = 1
    logqD = sapply(FUN = function(x) log(q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                           f0, f1, type, var_margy0, var_margy1, p)), D)
    for(i in 1:(N - 1)){
      for(j in (i + 1):N){
        pairwise_terms[counter] = k * logqD[i] + k * logqD[j] - k * log(sqrt((D[i] - D[j])^2))
        counter = counter + 1
      }
    }
    return(exp((1 / k) * logSumExp(pairwise_terms)))
  }
}


crit_fast = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / sqrt((D[i] - D[j])^2)
      counter = counter + 1
    }
  }
  return(max(pairwise_PEs))
}





evalD = function(D, N, p) det(crossprod(D))

evalDe = function(D, N, p) det(crossprod(D))^(1/p) / N

evalA = function(D, N, p) {
  Mi = solve(crossprod(D) / N)
  return(sum(diag(Mi))/p)
}

evalI = function(X, D, N, p) {
  Mi = solve(crossprod(D) / N)
  XMiX = rep(NA, dim(X)[1])
  for(i in 1:dim(X)[1]){
    XMiX[i] = X[i, ] %*% Mi %*% X[i, ]
  }
  return(mean(XMiX))
}

#evalGe = function(X, D, N, p){
#  Mi = solve(crossprod(D) / N)
#  XMiX = rep(NA, dim(X)[1])
#  for(i in 1:dim(X)[1]){
#    XMiX[i] = X[i, ] %*% Mi %*% X[i, ]
#  }
#  return(p / max(XMiX))
#}




