### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

## Like bayes_linear_regression.R, but I try to implement some of the techniques in
##  "Deterministic Sampling of Expensive Posteriors Using Minimum Energy Designs,"
##  Joseph 2018

library(transport)
library(mined)

### Helper  Functions ###

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2){
  return(sqrt(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
}

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0, beta0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1, beta1
  var = var_e + x^2 * var_mean # variance of marginal dist of y | H1, betai
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
# select the one with the smallest value of f_min to add to current design set D.
#  (Notice that this is different from the 2015 version, bc take MAX of sapply instead of SUM of sapply.)
f_min = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, 
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
    return(list("beta0" = beta0, "beta1" = beta1, "D" = D, "candidates" = C))
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
      
      f_min_candidates = sapply(C[[k]], function(x) f_min(x, D[1:(j - 1), k], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k] = C[[k]][chosen_cand]
      #Wass_D[j, k] = Wasserstein_distance(mean_beta0 * D[j, k], mean_beta1 * D[j, k], var_e + D[j, k]^2 * var_mean, var_e + D[j, k]^2 * var_mean)
    }
  }
  return(list("beta0" = beta0, "beta1" = beta1, "D" = D, "candidates" = C))
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 1 # same variance

n = 11 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5p = 103
K = 16
p = 1
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, xmin, xmax, K, p)

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model


# pictures

test_k = 1
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

#dev.copy(png,'efficient_bayes_lm_1.png')
#dev.off()

test_k = 4
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

#dev.copy(png,'efficient_bayes_lm_4.png')
#dev.off()

test_k = 8
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

#dev.copy(png,'efficient_bayes_lm_8.png')
#dev.off()

test_k = 16
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

#dev.copy(png,'efficient_bayes_lm_16.png')
#dev.off()


length(X_test$candidates[[3]]) # why is it 132?



