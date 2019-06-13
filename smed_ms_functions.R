###############################################################################################
# --- Functions Used in both One-At-a-Time and Fast Algorithms for Linear Model Selection --- #
###############################################################################################

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

# Calculate Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly = function(x, var_e, var_mean) var_e + x^2 * var_mean

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly(x, var_e, var_mean) # variance of marginal dist of y | H0 or H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var) # What about different vars? ************************************
  return(1.0 / Wass_dist^(1/2))  # ************** change from 1/2 to 1/(2p) at some point!
  # What does this imply the distribution of points is aymptotically? ************************************
}


############################################################################
# --- Functions for One-At-a-Time Algorithm for Linear Model Selection --- #
############################################################################

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
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = 10^5, 
                   k = 4, xmin = 0, xmax = 1, plotD = FALSE, genCandidates = 1, initialpt = 1){
  #  f0, mean_beta0 : regression line and mean slope of null model
  #  f1, mean_beta1 : regression line and mean slope of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on inputs
  # genCandidates :
  #   if 1, candidates are sequence from xmin to xmax with length numCandidates
  #   if 2, candidates are uniformly generated from xmin to xmax
  # initialpt : 
  #   how to get first design point in D
  #   if 1, minimize q
  #   if 2, from candidates (for more expensive functions?)minimize q
  
  # Draw a slope for each model
  beta0 = rnorm(n = 1, mean = mean_beta0, sd = sqrt(var_mean))
  beta1 = rnorm(n = 1, mean = mean_beta1, sd = sqrt(var_mean))
  
  # Create linear model
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  # Why are these needed? ************************************
  
  
  # -- Initialize 1st Design Point in D -- #
  # Joseph et al. equation (9) says to minimize q(x) to pick the first x ************************************
  D = rep(NA, N)
  if(initialpt == 2){
    xinitind = which.max(abs(f0(candidates) - f1(candidates)))
    D[1] = candidates[xinitind] # x1, first element of set of design points, D
    # candidates that are leftover (options to choose from for next 2:N design points)
    # candidates = candidates[-xinitind] # candidate set, for choosing next design point x_{n+1}
  } else{
    D[1] = optimize(function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  }
  
  # Plot density and (highest lik) points
  if(plotD == TRUE){
    curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
    curve(f1, col = 1, add = TRUE)
    text(D, f0(D), 1, col=4)
    text(D, f1(D), 1, col=4)
    points(D, 0, col=2)
  }

  # Sequentially pick rest
  # Now you probably do need the candidates because the problem is not convex ************************************
  # How long does it take? Can we use 10^5 candidates?  ************************************
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min(x, D[1:(i - 1)], k, mean_beta0, mean_beta1, var_e, var_mean))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew # add the new point to the set
    #candidates = candidates[-f_opt]
    if(plotD == TRUE){
      text(xnew, f0(xnew), i, col = 4)
      text(xnew, f1(xnew), i, col = 4)
      points(xnew, 0, col = 2)
    }
  }
  return(D)
}





###################################################################
# --- Functions for Fast Algorithm for Linear Model Selection --- #
###################################################################

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D.
#  (Notice that this is different from the 2015 version, bc take MAX of sapply instead of SUM of sapply.)
f_min_fast = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

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



### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 

# --- Functions for Fast Algorithm for Linear Model Selection --- #

SMED_ms_fast = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = NULL, 
                        xmin = 0, xmax = 1, K, p = 2, genCandidates = 1, initialpt = 1){
  #  N : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  numCandidates : number of candidates to generate at each step, if NULL, set to N by default.
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  # genCandidates :
  #   if 1, generate candidates in each neighborhood of each point using sequence with numCandidates elements
  #   if 2, generate uniformly in each neighborhood of each point
  #   if 3, generate using Lattice function, like they did...
  # initialpt : 
  #   how to get first design point in D
  #   if 1, minimize q
  #   if 2, from candidates
  
  if(is.null(numCandidates)) numCandidates = N
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
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
  
  optimal_q = optimize(function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  
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
      w_evals = sapply(C[[1]], FUN = function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                      var_marginaly(x, var_e, var_mean), 
                                                                      var_marginaly(x, var_e, var_mean)))
      # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
      xinitind = which.max(w_evals) # Here you could just use optimize ************************************
      #############   
      ############
      ############
      ############
      ############
      ############
      ############
      
      D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
      
      #Wass_D[1, k] = ... # a next step would be to matricize/preallocate these values for faster computing! ********************
      
    } else{
      D[1, k + 1] = optimal_q
    }
    
    
    # for j = 2:N
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # Shouldn't this be -R_jk ... and 0 not 1? ************************************fixed
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        # Did you check their code to verify this is how to construct the candidate set? ************************************
        # (the paper doesn't seem very clear) ************************************
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = p) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        
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






####################################################################################
####################################################################################
####################################################################################
####################################################################################
### Mixed Algorithms ###############################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################





# --- One-at-a-Time Criterion with Fast Method --- #

# have points converge to those which minimize the one-at-a-time algorithm's 
#  criterion (summation raised to k = 4 power) over K designs


f_min_1attimecrit_fast = function(candidate_jk, D_k, gamma_k, k_power, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^(gamma_k * k) * 
    sum(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))^k_power))
}

SMED_ms_1attimecrit_fast = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = NULL, 
                                    xmin = 0, xmax = 1, K, k_power = 4, p = 2, genCandidates = 1, initialpt = 1){
  #  N : number of design points to select (for set of design points, D_k each k = 1, ..., K)
  #  numCandidates : number of candidates to generate at each step, if NULL, set to N by default.
  #  K : number of designs to make (iteratively)
  #  xmin, xmax : limits on inputs
  # genCandidates :
  #   if 1, generate candidates in each neighborhood of each point using sequence with numCandidates elements
  #   if 2, generate uniformly in each neighborhood of each point
  #   if 3, generate using Lattice function, like they did...
  # initialpt : 
  #   how to get first design point in D
  #   if 1, minimize q
  #   if 2, from candidates
  
  if(is.null(numCandidates)) numCandidates = N
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = rep(NA, numCandidates)
  if(genCandidates == 1) C1 = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) C1 = sort(runif(numCandidates, xmin, xmax))
  if(genCandidates == 3) C1 = mined::Lattice(numCandidates, p = 1)
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
  
  optimal_q = optimize(function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  
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
      if(genCandidates == 3) tildeD1_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
      # save the candidates to be used in future designs
      C[[1]] = c(C[[1]], tildeD1_kplus)
      # criterion to choose first candidate from candidate set: 
      # the point at which f1 and f2 are most different
      w_evals = sapply(C[[1]], FUN = function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                      var_marginaly(x, var_e, var_mean), 
                                                                      var_marginaly(x, var_e, var_mean)))
      # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
      xinitind = which.max(w_evals) # Here you could just use optimize ************************************
      #############   
      ############
      ############
      ############
      ############
      ############
      ############
      
      D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
      
      #Wass_D[1, k] = ... # a next step would be to matricize/preallocate these values for faster computing! ********************
      
    } else{
      D[1, k + 1] = optimal_q
    }
    
    
    # for j = 2:N
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # Shouldn't this be -R_jk ... and 0 not 1? ************************************fixed
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        # Did you check their code to verify this is how to construct the candidate set? ************************************
        # (the paper doesn't seem very clear) ************************************
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_1attimecrit_fast(x, D[1:(j - 1), k + 1], gammas[k], k_power, mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_1attimecrit_fast(x, D[1:(j - 1), k + 1], gammas[k], k_power, mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}






f_min_fastcrit_1attime = function(candidate, D, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate, mean_beta0, mean_beta1, var_e, var_mean) * 
    max(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean) / abs(x_i - candidate))^k))
  # We might do a multivariate example at some point. 
}

# --- Fast Criterion with One-at-a-Time Method --- #

# sequentially pick each point (using a greedy algorithm) but by minimizing the fast
#  algorithm's criterion instead (max instead of summation)
SMED_ms_fastcrit_1attime = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = 10^5, 
                   xmin = 0, xmax = 1, plotD = FALSE, genCandidates = 1, initialpt = 1){
  #  f0, mean_beta0 : regression line and mean slope of null model
  #  f1, mean_beta1 : regression line and mean slope of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on inputs
  # genCandidates :
  #   if 1, candidates are sequence from xmin to xmax with length numCandidates
  #   if 2, candidates are uniformly generated from xmin to xmax
  # initialpt : 
  #   how to get first design point in D
  #   if 1, minimize q
  #   if 2, from candidates (for more expensive functions?)minimize q
  
  # Draw a slope for each model
  beta0 = rnorm(n = 1, mean = mean_beta0, sd = sqrt(var_mean))
  beta1 = rnorm(n = 1, mean = mean_beta1, sd = sqrt(var_mean))
  
  # Create linear model
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  # Why are these needed? ************************************
  
  
  # -- Initialize 1st Design Point in D -- #
  # Joseph et al. equation (9) says to minimize q(x) to pick the first x ************************************
  D = rep(NA, N)
  if(initialpt == 2){
    xinitind = which.max(abs(f0(candidates) - f1(candidates)))
    D[1] = candidates[xinitind] # x1, first element of set of design points, D
    # candidates that are leftover (options to choose from for next 2:N design points)
    # candidates = candidates[-xinitind] # candidate set, for choosing next design point x_{n+1}
  } else{
    D[1] = optimize(function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  }
  
  # Plot density and (highest lik) points
  if(plotD == TRUE){
    curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
    curve(f1, col = 1, add = TRUE)
    text(D, f0(D), 1, col=4)
    text(D, f1(D), 1, col=4)
    points(D, 0, col=2)
  }
  
  # Sequentially pick rest
  # Now you probably do need the candidates because the problem is not convex ************************************
  # How long does it take? Can we use 10^5 candidates?  ************************************
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_fastcrit_1attime(x, D[1:(i - 1)], mean_beta0, mean_beta1, var_e, var_mean))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew # add the new point to the set
    #candidates = candidates[-f_opt]
    if(plotD == TRUE){
      text(xnew, f0(xnew), i, col = 4)
      text(xnew, f1(xnew), i, col = 4)
      points(xnew, 0, col = 2)
    }
  }
  return(D)
}

































####################################################################################
####################################################################################
####################################################################################
####################################################################################
### Unknown Intercept ##############################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

var_marginaly2 = function(x, var_e, var_mean) var_e + (x^2 + 1) * var_mean

q2 = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly2(x, var_e, var_mean) # variance of marginal dist of y | H0 or H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}

### one at time ###

f_min2 = function(candidate, D, k, mean_beta0, mean_beta1, var_e, var_mean){
  q2(candidate, mean_beta0, mean_beta1, var_e, var_mean)^k * 
    sum(sapply(D, function(x_i) (q2(x_i, mean_beta0, mean_beta1, var_e, var_mean) / abs(x_i - candidate))^k))
}

SMED_ms2 = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = 10^5, 
                   k = 4, xmin = 0, xmax = 1, plotD = FALSE, genCandidates = 1, initialpt = 1){
  
  # Draw a slope for each model
  beta0 = rnorm(n = 1, mean = mean_beta0, sd = sqrt(var_mean))
  beta1 = rnorm(n = 1, mean = mean_beta1, sd = sqrt(var_mean))
  
  # Create linear model
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Generate Candidate Points -- #
  if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  
  # -- Initialize 1st Design Point in D -- #
  D = rep(NA, N)
  if(initialpt == 2){
    xinitind = which.max(abs(f0(candidates) - f1(candidates)))
    D[1] = candidates[xinitind]
  } else{
    D[1] = optimize(function(x) q2(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  }
  
  # Plot density and (highest lik) points
  if(plotD == TRUE){
    curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
    curve(f1, col = 1, add = TRUE)
    text(D, f0(D), 1, col=4)
    text(D, f1(D), 1, col=4)
    points(D, 0, col=2)
  }
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min2(x, D[1:(i - 1)], k, mean_beta0, mean_beta1, var_e, var_mean))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
    if(plotD == TRUE){
      text(xnew, f0(xnew), i, col = 4)
      text(xnew, f1(xnew), i, col = 4)
      points(xnew, 0, col = 2)
    }
  }
  return(D)
}


### fast ###

f_min_fast2 = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q2(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q2(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}


SMED_ms_fast2 = function(mean_beta0, mean_beta1, var_e, var_mean, N = 11, numCandidates = NULL, 
                         xmin = 0, xmax = 1, K, p = 2, genCandidates = 1, initialpt = 1){
  
  if(is.null(numCandidates)) numCandidates = N
  
  ## -- Create linear models -- #
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = rep(NA, numCandidates)
  if(genCandidates == 1) C1 = seq(from = xmin, to = xmax, length.out = numCandidates)
  if(genCandidates == 2) C1 = sort(runif(numCandidates, xmin, xmax))
  if(genCandidates == 3) C1 = mined::Lattice(numCandidates, p = 1)
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
  
  optimal_q = optimize(function(x) q2(x, mean_beta0, mean_beta1, var_e, var_mean), interval = c(xmin, xmax))$minimum
  
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
      if(genCandidates == 3) tildeD1_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
      # save the candidates to be used in future designs
      C[[1]] = c(C[[1]], tildeD1_kplus)
      # criterion to choose first candidate from candidate set: 
      # the point at which f1 and f2 are most different
      w_evals = sapply(C[[1]], FUN = function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                      var_marginaly2(x, var_e, var_mean), 
                                                                      var_marginaly2(x, var_e, var_mean)))
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
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # Shouldn't this be -R_jk ... and 0 not 1? ************************************fixed
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        # Did you check their code to verify this is how to construct the candidate set? ************************************
        # (the paper doesn't seem very clear) ************************************
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast2(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus = rep(NA, numCandidates)
        if(genCandidates == 1) tildeDj_kplus = seq(from = lower, to = upper, length.out = numCandidates)
        if(genCandidates == 2) tildeDj_kplus = runif(numCandidates, lower, upper)
        if(genCandidates == 3) tildeDj_kplus = mined::Lattice(numCandidates, p = 1) * (upper - lower) + lower
        C[[j]] = c(C[[j]], tildeDj_kplus) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast2(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
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
### Functions for Evaluating Design ################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################


### --- Variance on Slope --- ###

varslope_fixedbeta = function(D, var_e){
  xbar = mean(D)
  wi_denom = sum((D - xbar)^2)
  wi = sapply(D, FUN = function(x) (x - xbar)) / wi_denom
  return(var_e * sum(wi^2))
}

varslope = function(D, var_e, var_mean){
  xbar = mean(D)
  wi_denom = sum((D - xbar)^2)
  wi = sapply(D, FUN = function(x) (x - xbar)) / wi_denom
  var_marginaly_vec = sapply(D, FUN = function(x) var_marginaly(x, var_e, var_mean))
  return(sum(wi^2 * var_marginaly_vec))
}

varslope2 = function(D, var_e, var_mean){ # same result as varslope, just written s.t. first term is varslope_fixedbeta result
  xbar = mean(D)
  wi_denom = sum((D - xbar)^2)
  wi = sapply(D, FUN = function(x) (x - xbar)) / wi_denom
  first_term = var_e / wi_denom # same as varslope_fixedbeta
  second_term = var_mean * sum(D^2 * sapply(D, FUN = function(x) (x - xbar))^2) / wi_denom^2
  return(first_term + second_term)
}


### --- Compute Criterion --- ###

totalPE = function(D, N, mean_beta0, mean_beta1, var_e, var_mean){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / (D[i] - D[j])^2
      counter = counter + 1
    }
  }
  return(sum(pairwise_PEs))
}


crit_1atatime = function(D, N, k, mean_beta0, mean_beta1, var_e, var_mean){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = (qD[i] * qD[j] / (D[i] - D[j])^2)^k
      counter = counter + 1
    }
  }
  return((sum(pairwise_PEs))^(1/k))
}

 
crit_fast = function(D, N, mean_beta0, mean_beta1, var_e, var_mean){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / (D[i] - D[j])^2
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

evalGe = function(X, D, N, p){
  Mi = solve(crossprod(D) / N)
  XMiX = rep(NA, dim(X)[1])
  for(i in 1:dim(X)[1]){
    XMiX[i] = X[i, ] %*% Mi %*% X[i, ]
  }
  return(p / max(XMiX))
}


