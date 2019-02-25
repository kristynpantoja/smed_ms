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
isprime <- function(num) {
  if (num == 2) {
    TRUE
  } else if (any(num %% 2:(num-1) == 0)) {
    FALSE
  } else { 
    TRUE
  }
}

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    sum(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, numCandidates = 1000, 
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
  # check that n == numCandidates
  if(n != numCandidates) warning("number of samples (n) does not equal number of candidates (numCandidates).")
  # check that n >= 3
  if(n < 3) stop("not enough samples - need at least 3.")
  # check that numCandidates is a prime number
  if(!isprime(numCandidates)) stop("numCandidates is not prime!")
  candidates = mined::Lattice(numCandidates, p = 1)
  
  # -- Initialize Design D = {x_1, ..., x_n} -- #
  D = sort(candidates)
  # calculat logf for candidates
  logf0 = sapply(D, logf0)
  logf1 = sapply(D, logf1)
  
  # Plot initial design
  curve(f0, from = xmin, to = xmax)
  curve(f1, col = 2, add = TRUE)
  points(D, rep(0, numCandidates), col=3)
  
  # -- If K > 1, choose new design -- #
  if(K == 1) return(D)
  D = matrix(rep(D, K), nrow = K, ncol = n, byrow = T)
  gammas = (c(1:K) - 1) / K
  for(k in 2:K){
    for(j in 2:(n - 1)){
      # get candidates in neighborhood L_jk = (lower, upper)
      R_jk = which.min(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k])
      lower = D[j, k] - R_jk
      upper = D[j, k] + R_jk
      candidates_jk = Lattice(n, p = 1) * (upper - lower) + lower
      # calculate log f-hat of candidates_jk (!!!!!! for now, f-hat = f)
      cand_logf0 = sapply(candidates_jk, logf0)
      cand_logf1 = sapply(candidates_jk, logf1)
      # pick from candidates
      # calculate criterion for each point
      f_min_candidates = sapply(candidates_jk, 
                                function(x) f_min(x, D[1:(j - 1), k], gammas[k], 
                                                  mean_beta0, mean_beta1, 
                                                  var_e, var_mean))
      D[j, k] = which.min(f_min_candidates)
    }
  }
  return(D)
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.00001
var_e = 1 # same variance

n = 7 # in paper, n = numCandidates
numCandidates = 7 # largest prime number less than 100 + 5p = 103
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)
