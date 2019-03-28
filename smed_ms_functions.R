# --- Functions Used in both One-At-a-Time and Fast Algorithms for Linear Model Selection --- #

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2){
  return(sqrt(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
}

var_marginaly = function(x, var_e, var_mean) var_e + x^2 * var_mean

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly(x, var_e, var_mean) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}



# --- Functions for One-At-a-Time Algorithm for Linear Model Selection --- #

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min = function(candidate, D, k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate, mean_beta0, mean_beta1, var_e, var_mean)^k * 
    sum(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean) / abs(x_i - candidate))^k))
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
  candidates = runif(numCandidates, xmin, xmax)
  
  # -- Initialize 1st Design Point in D -- #
  # get the point at which f1 and f2 are most different
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

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D.
#  (Notice that this is different from the 2015 version, bc take MAX of sapply instead of SUM of sapply.)
f_min_fast = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms_fast = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
                        xmin = 0, xmax = 1, K, p){
  #  n : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  numCandidatesinit = n
  C1 = mined::Lattice(numCandidatesinit, p = 1)
  D1 = sort(C1)
  C = C1
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K = 1, return the design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = n, ncol = K)
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = (c(1:K) - 1) / (K - 1)
  # save candidates for each K
  C <- list()
  for (k in 1:K){
    C[[k]] = C1
  }
  
  for(k in 2:K){
    # for j = 1
    # get candidates in L_1k
    R1k = D[2, k] - D[1, k] # radius of L1k
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ################################
    L1k_upper = min(D[1, k] + R1k, 1)
    # candidates from space-filling design, tildeD_1k
    tildeD1k = Lattice(numCandidatesinit, p = 1) * (L1k_upper - L1k_lower) + L1k_lower
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    # criterion to choose candidate from candidate set - for now, choose like in older paper:
    # get the point at which f1 and f2 are most different
    C[[k]] = c(C[[k]], tildeD1k)
    xinitind = which.max(abs(f0(C[[k]]) - f1(C[[k]])))
    D[1, k] = C[[k]][xinitind] # x1, first element of set of design points, D
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:n){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == n){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = D[j, k] - R_jk
        upper = D[j, k] + R_jk
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      }
      
      # calculate Wassersteins for each candidate
      #Wass_cand = sapply(C[[k]], function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
      # pick from candidates:
      # first, calculate criterion for each candidate
      
      f_min_candidates = sapply(C[[k]], function(x) f_min_fast(x, D[1:(j - 1), k], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = C[[k]][chosen_cand]
      #Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each



# --- Functions for Fast Algorithm for Linear Model Selection,
#       to Check Optimization and See Why Results Aren't As Expected,
#       i.e. Not Like One-At-a-Time Algorithm for LMS --- #


# This is somewhat a hybrid between f_min (with sum & k_power) and f_min_fast (with max & gamma_k)
#  in that now, instead of max & gamma_k as in f_min_fast, we do
#  sum with both k_power and gamma_k.
f_min_fast2 = function(candidate_jk, D_k, gamma_k, k_power, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^(k_power * gamma_k) * 
    sum(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k 
                                   / abs(x_i - candidate_jk))^k_power))
}

SMED_ms_fast2 = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
                        xmin = 0, xmax = 1, K, p, k_power){
  #  n : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  numCandidatesinit = n
  C1 = mined::Lattice(numCandidatesinit, p = 1)
  D1 = sort(C1)
  C = C1
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K = 1, return the design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = n, ncol = K)
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = (c(1:K) - 1) / (K - 1)
  # save candidates for each K
  C <- list()
  for (k in 1:K){
    C[[k]] = C1
  }
  
  for(k in 2:K){
    # for j = 1
    # get candidates in L_1k
    R1k = D[2, k] - D[1, k] # radius of L1k
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ################################
    L1k_upper = min(D[1, k] + R1k, 1)
    # candidates from space-filling design, tildeD_1k
    tildeD1k = Lattice(numCandidatesinit, p = 1) * (L1k_upper - L1k_lower) + L1k_lower
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    # criterion to choose candidate from candidate set - for now, choose like in older paper:
    # get the point at which f1 and f2 are most different
    C[[k]] = c(C[[k]], tildeD1k)
    xinitind = which.max(abs(f0(C[[k]]) - f1(C[[k]])))
    D[1, k] = C[[k]][xinitind] # x1, first element of set of design points, D
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:n){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == n){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = D[j, k] - R_jk
        upper = D[j, k] + R_jk
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      }
      
      # calculate Wassersteins for each candidate
      #Wass_cand = sapply(C[[k]], function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
      # pick from candidates:
      # first, calculate criterion for each candidate
      
      f_min_candidates = sapply(C[[k]], function(x) f_min_fast2(x, D[1:(j - 1), k], gammas[k], k_power, mean_beta0, mean_beta1, var_e, var_mean))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = C[[k]][chosen_cand]
      #Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}




# --- Functions for Fast Algorithm for Linear Model Selection,
#      Without min-ing and max-ing lower and upper bounds of Ljk --- #

SMED_ms_fast3 = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
                        xmin = 0, xmax = 1, K, p){
  #  n : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  numCandidatesinit = n
  C1 = mined::Lattice(numCandidatesinit, p = 1)
  D1 = sort(C1)
  C = C1
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K = 1, return the design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = n, ncol = K)
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = (c(1:K) - 1) / (K - 1)
  # save candidates for each K
  C <- list()
  for (k in 1:K){
    C[[k]] = C1
  }
  
  for(k in 2:K){
    # for j = 1
    # get candidates in L_1k
    R1k = D[2, k] - D[1, k] # radius of L1k
    L1k_lower = D[1, k] - R1k
    L1k_upper = D[1, k] + R1k
    # candidates from space-filling design, tildeD_1k
    tildeD1k = Lattice(numCandidatesinit, p = 1) * (L1k_upper - L1k_lower) + L1k_lower
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    # criterion to choose candidate from candidate set - for now, choose like in older paper:
    # get the point at which f1 and f2 are most different
    C[[k]] = c(C[[k]], tildeD1k)
    xinitind = which.max(abs(f0(C[[k]]) - f1(C[[k]])))
    D[1, k] = C[[k]][xinitind] # x1, first element of set of design points, D
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:n){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == n){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = D[j, k] - R_jk
        upper = D[j, k] + R_jk
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      }
      
      # calculate Wassersteins for each candidate
      #Wass_cand = sapply(C[[k]], function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
      # pick from candidates:
      # first, calculate criterion for each candidate
      
      f_min_candidates = sapply(C[[k]], function(x) f_min_fast(x, D[1:(j - 1), k], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = C[[k]][chosen_cand]
      #Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}







# --- Functions for Fast Algorithm for Linear Model Selection,
#       to Check Optimization --- #


# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D.
#  (Notice that this is different from the 2015 version, bc take MAX of sapply instead of SUM of sapply.)
f_min_fast4 = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms_fast4 = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
                        xmin = 0, xmax = 1, K, p){
  #  n : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  # check that n is the largest prime number less than 100 + 5p.... ###################################
  # or at least check that it's prime!
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  numCandidatesinit = n
  C1 = mined::Lattice(numCandidatesinit, p = 1)
  D1 = sort(C1)
  C = C1
  
  ## calculate Wassersteins for each design point
  #Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K = 1, return the space-filling design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = n, ncol = K)
  #Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = (c(1:K) - 1) / (K - 1)
  # save candidates for each K
  C <- list()
  for (k in 1:K){
    C[[k]] = C1 # = tilde_D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # Get candidates in neighborhood L1k = (lower, upper):
    # 1. find point with smallest distance from 1st point, D[1, k]:
    differences = sapply(D[ , k], FUN = function(x) abs(D[1, k] - x))
    differences[1] = NA # so that 1st point isn't selected
    # 2. calculate radius of L1k, R1k = distance between D[1, k]
    R1k = differences(which.min(differences)) # radius of L1k
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? ####################
    L1k_upper = min(D[1, k] + R1k, 1)
    # candidates from space-filling design, tildeD_1k
    tildeD1k = Lattice(numCandidatesinit, p = 1) * (L1k_upper - L1k_lower) + L1k_lower
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    # criterion to choose candidate from candidate set - for now, choose like in older paper:
    # get the point at which f1 and f2 are most different
    C[[k]] = c(C[[k]], tildeD1k)
    xinitind = which.max(abs(f0(C[[k]]) - f1(C[[k]])))
    D[1, k] = C[[k]][xinitind] # x1, first element of set of design points, D
    #Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:n){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == n){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = D[j, k] - R_jk
        upper = D[j, k] + R_jk
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        C[[k]] = c(C[[k]], tildeDjk)
      }
      
      # calculate Wassersteins for each candidate
      #Wass_cand = sapply(C[[k]], function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, var_e + x^2 * var_mean, var_e + x^2 * var_mean))
      # pick from candidates:
      # first, calculate criterion for each candidate
      
      f_min_candidates = sapply(C[[k]], function(x) f_min_fast(x, D[1:(j - 1), k], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = C[[k]][chosen_cand]
      #Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}






