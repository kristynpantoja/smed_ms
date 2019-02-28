### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

## Like bayes_linear_regression.R, but I try to implement some of the techniques in
##  "Deterministic Sampling of Expensive Posteriors Using Minimum Energy Designs,"
##  their 2018 paper.

library(transport)
library(mined)

### Helper  Functions ###

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var_1, var_2){
  return(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2)
}

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_e + x^2 * var_mean # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
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
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
                   xmin = 0, xmax = 1, K, p){
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
  f0 = function(x) beta0 * x # null regression model
  f1 = function(x) beta1 * x # alternative regression model
  # log models - didn't use!
  logf0 = function(x) log(beta0) + log(x) # null regression model
  logf1 = function(x) log(beta1) + log(x) # alternative regression model
  
  # -- Generate Candidate Points -- #
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  
  numCandidatesinit = 7 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!# largest prime number less than 100 + 5p = 103 but MAKE IT A FUNCTION.
  candidates = mined::Lattice(numCandidatesinit, p = 1)
  
  # -- Initialize Design D = {x_1, ..., x_n} -- #
  D_init = sort(candidates)
  # calculate logf for candidates
  logf0 = sapply(D_init, logf0)
  logf1 = sapply(D_init, logf1)
  # calculate Wassersteins for each design point
  Wasser_init = sapply(D_init, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                var_e + x^2 * var_mean, var_e + x^2 * var_mean))
  
  # -- If K > 1, choose new design -- #
  if(K == 1) return(D)
  D = matrix(rep(D_init, K), nrow = n, ncol = K)
  Wass_D = matrix(rep(Wasser_init, K), nrow = n, ncol = K)
  gammas = (c(1:K) - 1) / (K - 1)
  for(k in 2:K){
    # for j = 1
    # get candidates in R_1k
    R_1k = D[2, k] - D[1, k]
    lower = max(D[1, k] - R_1k, 0) # is this necessary, to max with 0? ################################
    upper = min(D[1, k] + R_1k, 1)
    candidates_1k = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
    # save the candidates to be used in future designs
    candidates_1k = c(candidates_1k, candidates)
    candidates = candidates_1k
    # criterion to choose candidate from candidate set - for now, choose like in older paper:
    # get the point at which f1 and f2 are most different
    xinitind = which.max(abs(f0(candidates_1k) - f1(candidates_1k)))
    D[1, k] = candidates_1k[xinitind] # x1, first element of set of design points, D
    Wass_D[1, k] = Wasserstein_distance(mean_beta0 * xinitind, mean_beta1 * xinitind, 
                                        var_e + xinitind^2 * var_mean, var_e + xinitind^2 * var_mean)
    # at some point... see their code for initializing kth design. ##################################
    
    # for j = 2:n
    for(j in 2:n){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == n){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = D[j, k] - R_jk
        upper = D[j, k] + R_jk
        candidates_jk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        candidates_jk = c(candidates_jk, candidates)
        candidates = candidates_jk
      } else{
        R_jk = which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) # is this necessary, to max with 0? ################################
        upper = min(D[j, k] + R_jk, 1)
        candidates_jk = Lattice(numCandidatesinit, p = 1) * (upper - lower) + lower
        candidates_jk = c(candidates_jk, candidates)
        candidates = candidates_jk
      }
      numCandidates = length(candidates_jk)
      # calculate log f-hat of candidates_jk (!!!!!! for now, f-hat = f)
      # cand_logf0 = sapply(candidates_jk, logf0)
      # cand_logf1 = sapply(candidates_jk, logf1)
      
      # calculate Wassersteins for each candidate
      Wass_cand = sapply(candidates_jk, function(x) Wasserstein_distance(mean_beta0 * x, mean_beta1 * x, 
                                                                         var_e + x^2 * var_mean, var_e + x^2 * var_mean))
      # pick from candidates:
      # first, calculate criterion for each candidate
      f_min_candidates = rep(NA, numCandidates)
      for(m in 1:numCandidates){
        f_min_candidates[m] = Wass_cand[m]^(gammas[k] / 2) * sum((Wass_D[1:(j - 1), k]^(gammas[k] / 2)) * abs(D[1:(j - 1), k] - candidates_jk[m])^(2))
      }
      
      
      #f_opt = which.min(f_min(candidates, D, k, mean_beta0, mean_beta1, var_e, var_mean))
      f_min_candidates = sapply(candidates, function(x) f_min(x, D[1:(j - 1), k], 4, mean_beta0, mean_beta1, var_e, var_mean))
      
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = candidates_jk[chosen_cand]
      Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], 
                                          var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = beta0, "beta1" = beta1, "D" = D, "candidates" = candidates))
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 1 # same variance

n = 7 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5p = 103
K = 4
p = 1
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, xmin, xmax, K, p)

f0 = function(x) X_test$beta0 * x # null regression model
f1 = function(x) X_test$beta1 * x # alternative regression model

test_k = 2
curve(f0, from = xmin, to = xmax)
curve(f1, col = 2, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

#sort(X_test$candidates)
length(X_test$candidates)



