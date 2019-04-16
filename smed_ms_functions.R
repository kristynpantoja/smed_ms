# --- Functions Used in both One-At-a-Time and Fast Algorithms for Linear Model Selection --- #

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2, type = 1){
  if(type == 1){ # univariate case
    #wass = sqrt((mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
    wass = sqrt((mu1 - mu2)^2 + (sqrt(var1) - sqrt(var2))^2)
  } else{
    if(type > 1){ # multivariate case 
      sqrt_var2 = sqrtm(var2)
      wass = sqrt(crossprod(mu1 - mu2) + sum(diag(var1 + var2 - 2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2))))
    } else{
      stop("invalid type")
    }
  }
  return(as.numeric(wass))
}

# Is this var(beta*x + epsilon) ? 
var_marginaly = function(x, var_e, var_mean) var_e + x^2 * var_mean

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly(x, var_e, var_mean) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var) # What about different vars?
  return(1.0 / Wass_dist^(1/2)) 
  # What does this imply the distribution of points is aymptotically?
}



# --- Functions for One-At-a-Time Algorithm for Linear Model Selection --- #

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min = function(candidate, D, k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate, mean_beta0, mean_beta1, var_e, var_mean)^k * 
    sum(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean) / abs(x_i - candidate))^k))
    # We might do a multivariate example at some point. 
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, numCandidates = 1000, 
                   k = 4, xmin = 0, xmax = 1){
  #  f0, mean_beta0 : regression line and mean slope of null model
  #  f1, mean_beta1 : regression line and mean slope of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on inputs
  
  # Draw a slope for each model
  beta0 = rnorm(n = 1, mean = mean_beta0, sd = sqrt(var_mean))
  beta1 = rnorm(n = 1, mean = mean_beta1, sd = sqrt(var_mean))
  
  # Create linear model
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # Calculate posterior means and variance for the mean beta_i given the data y - DON'T NEED THESE.
  # where m = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta_i) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # and s = ((1 / sigma_e^2) + (1 / var_mean^2))^(-1)
  # post_beta0 = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta0) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # post_beta1 = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta1) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # post_var = ((1 / sigma_e^2) + (1 / var_mean^2))^(-1)
  
  # -- Generate Candidate Points -- #
<<<<<<< HEAD
  candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
=======
  candidates = runif(numCandidates, xmin, xmax) # Why are these needed?
>>>>>>> 30a005c4592ad071d4a0e4b1741496b56580218f
  
  # -- Initialize 1st Design Point in D -- #
  # get the point at which f1 and f2 are most different
  # Why not just pick the optimal x rather than randomly generating
  # points and picking the best one? That is maximize (f0(x)-f1(x))^2
  # using optimize (or our knowledge about the problem)
  # Joseph et al. equation (9) says to minimize q(x) to pick the first x
  xinitind = which.max(abs(f0(candidates) - f1(candidates)))
  D = candidates[xinitind] # x1, first element of set of design points, D
  # also D is used to initialize greedy algorithm
  
  # candidates that are leftover (options to choose from for next 2:n design points)
  candidates = candidates[-xinitind] # candidate set, for choosing next design point x_{n+1}
  
  # Plot density and (highest lik) points
  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
  curve(f1, col = 1, add = TRUE)
  text(D, f0(D), 1, col=4)
  text(D, f1(D), 1, col=4)
  points(D, 0, col=2)
  
  # Sequentially pick rest
  # Now you probably do need the cnadidates because the problem is not convex
  # How long does it take? Can we use 10^5 cnadidates? 
  for(i in 2:n){
    # Find f_opt: minimum of f_min
    
    #f_opt = which.min(f_min(candidates, D, k, mean_beta0, mean_beta1, var_e, var_mean))
    f_min_candidates = sapply(candidates, function(x) f_min(x, D, k, mean_beta0, mean_beta1, var_e, var_mean))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D = c(D,xnew) # add the new point to the set
    #candidates = candidates[-f_opt]
    text(xnew, f0(xnew), i, col = 4)
    text(xnew, f1(xnew), i, col = 4)
    points(xnew, 0, col = 2)
  }
  return(D)
}



# --- Functions for Fast Algorithm for Linear Model Selection --- #

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D.
#  (Notice that this is different from the 2015 version, bc take MAX of sapply instead of SUM of sapply.)
f_min_fast = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk)^(2))))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 

# --- Functions for Fast Algorithm for Linear Model Selection --- #

SMED_ms_fast = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, 
                         xmin = 0, xmax = 1, K, p = 1){
  #  N : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = seq(from = xmin, to = xmax, length.out = N)
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
  gammas = c(1:K) / (K - 1) # Shouldn't this be 0:(K-1)/(k-1) ?
  # -- the final step should be gamma = 1 because then we optimize the correct criterion
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # Get candidates in neighborhood L1k = (lower, upper):
    # for j = 1
    # get candidates in L_1k
    R1k = min(abs(D[-1, k] - D[1, k])) # radius of L1k # can do it like this bc sorted D1, which was used to initialize D
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ################################
    L1k_upper = min(D[1, k] + R1k, 1) # HERE IT IS BECAUSE o/w GET NaNs in q evaluation!
    # candidates from space-filling design, tildeD1_kplus1
    # In general the number of local points does no need to be N
    # so I suggest introducing a N_L. You can set N_L = N for now
    # but we may decide to change it later.
    tildeD1_kplus1 = seq(from = L1k_lower, to = L1k_upper, length.out = N)
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    C[[1]] = c(C[[1]], tildeD1_kplus1)
    # criterion to choose first candidate from candidate set: 
    # the point at which f1 and f2 are most different
    w_evals = sapply(C[[1]], FUN = function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                    var_marginaly(x, var_e, var_mean), 
                                                                    var_marginaly(x, var_e, var_mean)))
    xinitind = which.max(w_evals) # Here you could just use optimize
    
    #xinitind = which.max(abs(f0(C[[1]]) - f1(C[[1]])))
    
    ### is this the same as the one with largest f? (here, Wasserstein^(1/2p))
    ### i.e. the one with smallest q?
    ### turns out, yes!:
    ###   a = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    ###   xinitind == which.min(a)
    # why are they equivalent, though? is wasserstein always ? 1?
    # must be, since q is basically 1 / wasserstein?
    #q_evals = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    #xinitind = which.min(q_evals)
    
    D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
    
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = min(D[j, k] + R_jk, 1) # Shouldn't this be -R_jk ... and 0 not 1?
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus1 = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDj_kplus1) # This is now C_j^{k+1}
        # Did you check their code to verify this is how to construct the candidate set?
        # (the paper doesn't seem very clear)
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDjk) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each



###########
# It would be better to input a function for defining the candidate set
# rather than writing a seperate SMED function for each case


# --- Fast Algorithm for Linear Model Selection
#      With Uniform instead of Seq --- #

SMED_ms_fast2 = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, 
                        xmin = 0, xmax = 1, K, p = 1){
  #  N : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = sort(runif(N, xmin, xmax))
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
  gammas = c(1:K) / (K - 1) # Shouldn't this be 0:(K-1)/(k-1) ?
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # Get candidates in neighborhood L1k = (lower, upper):
    # for j = 1
    # get candidates in L_1k
    R1k = min(abs(D[-1, k] - D[1, k])) # radius of L1k # can do it like this bc sorted D1, which was used to initialize D
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ################################
    L1k_upper = min(D[1, k] + R1k, 1) # HERE IT IS BECAUSE o/w GET NaNs in q evaluation!
    # candidates from space-filling design, tildeD1_kplus1
    tildeD1_kplus1 = runif(N, L1k_lower, L1k_upper)
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    C[[1]] = c(C[[1]], tildeD1_kplus1)
    # criterion to choose first candidate from candidate set: the point at which f1 and f2 are most different
    q_evals = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    xinitind = which.min(q_evals)
    #xinitind = which.max(abs(f0(C[[1]]) - f1(C[[1]])))
    
    ### is this the same as the one with largest f? (here, Wasserstein^(1/2p))
    ### i.e. the one with smallest q?
    ### turns out, yes!:
    ###   a = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    ###   xinitind == which.min(a)
    
    D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
    
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = min(D[j, k] + R_jk, 1) # Is this right?
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus1 = runif(N, lower, upper)
        C[[j]] = c(C[[j]], tildeDj_kplus1) # This is now C_j^{k+1}
        
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = runif(N, lower, upper)
        C[[j]] = c(C[[j]], tildeDjk) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}


# --- Fast Algorithm for Linear Model Selection
#      With Lattice instead of Seq --- #

# check if a number is prime, to use Lattice function in
isprime <- function(x) {
  if (x == 2) {
    TRUE
  } else if (any(x %% 2:(x - 1) == 0)) {
    FALSE
  } else { 
    TRUE
  }
}

SMED_ms_fast3 = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, 
                        xmin = 0, xmax = 1, K, p = 1){
  #  N : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = mined::Lattice(N, p = 1)
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
  gammas = c(1:K) / (K - 1) # Shouldn't this be 0:(K-1)/(k-1) ?
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # Get candidates in neighborhood L1k = (lower, upper):
    # for j = 1
    # get candidates in L_1k
    R1k = min(abs(D[-1, k] - D[1, k])) # radius of L1k # can do it like this bc sorted D1, which was used to initialize D
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ################################
    L1k_upper = min(D[1, k] + R1k, 1) # HERE IT IS BECAUSE o/w GET NaNs in q evaluation!
    # candidates from space-filling design, tildeD1_kplus1
    tildeD1_kplus1 = mined::Lattice(N, p = 1) * (L1k_upper - L1k_lower) + L1k_lower
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    C[[1]] = c(C[[1]], tildeD1_kplus1)
    # criterion to choose first candidate from candidate set: the point at which f1 and f2 are most different
    q_evals = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    xinitind = which.min(q_evals)
    #xinitind = which.max(abs(f0(C[[1]]) - f1(C[[1]])))
    
    ### is this the same as the one with largest f? (here, Wasserstein^(1/2p))
    ### i.e. the one with smallest q?
    ### turns out, yes!:
    ###   a = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    ###   xinitind == which.min(a)
    
    D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
    
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = min(D[j, k] + R_jk, 1) # Is this right?
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus1 = mined::Lattice(N, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus1) # This is now C_j^{k+1}
        
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = mined::Lattice(N, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDjk) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}















####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################






