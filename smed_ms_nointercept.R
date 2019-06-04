library(transport)
library(mined)
library(expm)
source("/Users/kristyn/Documents/research/smed_ms/smed_ms_functions.R")


var_marginaly2 = function(x, var_e, var_mean) var_e + (x^2 + 1) * var_mean

q2 = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0[1] + mean_beta0[2] * x # mean of marginal dist of y | H0
  mu2 = mean_beta1[1] + mean_beta1[2] * x # mean of marginal dist of y | H1
  var = var_marginaly2(x, var_e, var_mean) # variance of marginal dist of y | H0 or H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var, dim = 2)
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
  beta0 = c(rnorm(n = 1, mean = mean_beta0[1], sd = sqrt(var_mean[1, 1])), 
            rnorm(n = 1, mean = mean_beta0[2], sd = sqrt(var_mean[2, 2])))
  beta1 = c(rnorm(n = 1, mean = mean_beta1[1], sd = sqrt(var_mean[1, 1])), 
            rnorm(n = 1, mean = mean_beta1[2], sd = sqrt(var_mean[2, 2])))
  
  # Create linear model
  f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
  f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
  
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



















##############################################################################

mean_beta0 = c(1/4, 1) # slope of null model
mean_beta1 = c(0, 1 / 2) # slope of alternative model
var_mean = matrix(c(0.001, 0, 0, 0.001), 2, 2, byrow = T) # variance on beta
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model

N = 67
# for one-at-a-time algorithm:
numCandidates = 100 # suggested 10^5

###############################
# One-at-a-Time Algorithm k = 4
k = 4
X_1atatime = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)

par(mfrow = c(1, 1))
# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime[i], y[i], i, col=4)
}
points(X_1atatime, rep(0, N), col = 2, pch = "*")
lines(X_1atatime, y, col = 3)

hist(X_1atatime, freq = T, main = "1atT k=4")
###############################


###############################
# One-at-a-Time Algorithm k = 1
k = 1
X_1atatime_k1 = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)


###############################


# for fast algorithm:
K = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
X_fast = SMED_ms_fast2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, xmin = xmin, xmax = xmax, K = K, p = p)










mean_beta0 = c(1, -1) # slope of null model
mean_beta1 = c(0, 1) # slope of alternative model
var_mean = matrix(c(0.001, 0, 0, 0.001), 2, 2, byrow = T) # variance on beta
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model

N = 67
# for one-at-a-time algorithm:
numCandidates = 100 # suggested 10^5

###############################
# One-at-a-Time Algorithm k = 4
k = 4
X_1atatime2 = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)

par(mfrow = c(1, 1))
# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime2[i], y[i], i, col=4)
}
points(X_1atatime2, rep(0, N), col = 2, pch = "*")
lines(X_1atatime2, y, col = 3)

hist(X_1atatime2, freq = T, main = "1atT k=4")
###############################


###############################
# One-at-a-Time Algorithm k = 1
k = 1
X_1atatime_k1 = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)


###############################


# for fast algorithm:
K = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
X_fast = SMED_ms_fast2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, xmin = xmin, xmax = xmax, K = K, p = p)





























---
  
  ```{r unkint, cache = T}
mean_beta0 = c(1/4, 1) # slope of null model
mean_beta1 = c(0, 1 / 2) # slope of alternative model
var_mean = matrix(c(0.001, 0, 0, 0.001), 2, 2, byrow = T) # variance on beta
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model

N = 67
# for one-at-a-time algorithm:
numCandidates = 10^5 # suggested 10^5

# One-at-a-Time Algorithm k = 4
k = 4
X_1atatime = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)

# One-at-a-Time Algorithm k = 1
k = 1
X_1atatime_k1 = SMED_ms2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax)



# for fast algorithm:
K = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
X_fast = SMED_ms_fast2(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, var_e = var_e, var_mean = var_mean, N = N, xmin = xmin, xmax = xmax, K = K, p = p)
```

## One-at-a-Time, k = 4

```{r, fig.height = 6}
par(mfrow = c(1, 2))
# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime[i], y[i], i, col=4)
}
points(X_1atatime, rep(0, N), col = 2, pch = "*")
lines(X_1atatime, y, col = 3)
# histogram of design points
hist(X_1atatime, freq = T, main = "1atT k=4")
```

## One-at-a-Time, k = 1

```{r, fig.height = 6}
par(mfrow = c(1, 2))
# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=1")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime_k1[i], y[i], i, col=4)
}
points(X_1atatime_k1, rep(0, N), col = 2, pch = "*")
lines(X_1atatime_k1, y, col = 3)
# histogram of design points
hist(X_1atatime_k1, freq = T, main = "1atT k=1")
```

## Fast Algorithm (K = 20)

```{r, fig.height = 6}
par(mfrow = c(1, 2))
# design points' locations
X_fast_K = X_fast$D[ , K]
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = FALSE, main = "Fast K=20")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_fast$D[i ,K], y[i], i, col=4)
}
points(X_fast_K, rep(0, N), col = 2, pch = "*")
lines(X_fast_K, y, col = 3)
# histogram of design points
hist(X_fast_K, freq = T, main = "Fast K=20")
```










