###########################
# Making W a Distribution #
###########################

# --- Directory --- #
# Kristyn
# home = "/Users/kristyn/Documents/research/smed_ms/"
# Temp
home = "/Users/kristyn/Desktop/research/smed_ms/"

# --- Sources/Libraries --- #
med_fns = paste(home,"med_ms_functions.R",sep="")
source(med_fns)
library(expm)
library(matrixStats)

# --- Functions (for printing) --- #

MED_ms_print = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
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
    print(paste(i, " out of ", N))
  }
  return(D)
}

MED_ms_fast_print = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
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
  C <- matrix(rep(NA, numCandidates * K * N), numCandidates * K, N)
  for(j in 1:N){
    C[ , j] = c(C1, rep(NA, numCandidates * (K - 1)))
  }
  
  
  if(initialpt == 2){
    w_evals = sapply(C1, FUN = function(x) Wasserstein_distance(f0(x), f1(x), 
                                                                var_marginaly(x, var_e, var_mean0, type, var_margy0), 
                                                                var_marginaly(x, var_e, var_mean1, type, var_margy1)))
    # Joseph et al.2018, after equation (8), says to maximize f(x) to pick the first x (which for us is Wass dist)
    xinitind = which.max(w_evals)
    optimal_val = C1[xinitind]
  } else{
    optimal_val = optimize(function(x) q(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, 
                                         type, var_margy0, var_margy1, p), interval = c(xmin, xmax))$minimum
  }
  
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    D[1, k + 1] = optimal_val
    
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
      C[ (k * numCandidates + 1):((k + 1) * numCandidates), j] = tildeDj_kplus # This is now C_j^{k+1}
      
      f_min_candidates = sapply(C[ , j], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], 
                                                                mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                                                                f0, f1, type, var_margy0, var_margy1, p))
      #choose that which has largest evaluation of criterion
      chosen_cand = which.min(f_min_candidates)
      D[j, k + 1] = C[ , j][chosen_cand]
      print(paste(k, " out of ", (K - 1), " stages : ", j, " out of ", N, "points."))
    }
  }
  
  return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
}

# --- Design Generation --- #

# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 100
type = c(2, 2)
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 4
# design_1atT = MED_ms_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                      var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                      f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                      N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
#                      log_space = TRUE, genCandidates = 1, initialpt = 1)

# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 200
type = c(2, 2)
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 4
# design_1atT2 = MED_ms_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                            var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                            f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                            N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
#                            log_space = TRUE, genCandidates = 1, initialpt = 1)



# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 100
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# design_fast1 = MED_ms_fast_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                                 var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                                 f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                                 N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
#                                 genCandidates = 1, initialpt = 1)
# saveRDS(design_fast1, paste(home, "limit_distr/", "design_fast_simplelinear", "_N", N, "_S",S, ".rds",sep=""))

# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 200
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# design_fast2 = MED_ms_fast_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                                  var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                                  f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                                  N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
#                                  genCandidates = 1, initialpt = 1)
# saveRDS(design_fast2, paste(home, "limit_distr/", "design_fast_simplelinear", "_N", N, "_S",S, ".rds",sep=""))


# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 300
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# design_fast3 = MED_ms_fast_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                                  var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                                  f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                                  N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
#                                  genCandidates = 1, initialpt = 1)
# saveRDS(design_fast3, paste(home, "limit_distr/", "design_fast_simplelinear", "_N", N, "_S",S, ".rds",sep=""))




# charge function evaluated at x
# q = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
#   if(length(type) != 2) stop("type should be vector with length == 2")
#   mu1 = f0(x) # mean of marginal dist of y | H0
#   mu2 = f1(x) # mean of marginal dist of y | H1
#   var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
#   var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
#   Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
#   return(1.0 / (Wass_dist^(1/(2*p))))
# }
# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 300
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# design_fast3.2 = MED_ms_fast_print(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1,
#                                  var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e,
#                                  f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL,
#                                  N = N, K = S, p = p, xmin = xmin, xmax = xmax,
#                                  genCandidates = 1, initialpt = 1)
# saveRDS(design_fast3.2, paste(home, "icloud/limit_distr/", "design_fast_simplelinear", "_N", N, "_S",S, "rho1",".rds",sep=""))


# --- Design Load --- #

design_fast1 = readRDS(paste(home, "icloud/limit_distr/design_fast_simplelinear_N100_S5.rds", sep = ""))
design_fast2 = readRDS(paste(home, "icloud/limit_distr/design_fast_simplelinear_N200_S5.rds", sep = ""))
design_fast3 = readRDS(paste(home, "icloud/limit_distr/design_fast_simplelinear_N300_S5.rds", sep = ""))

# --- Distribution of W --- #

# Model Parameters / Settings
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 300
type = c(2, 2)
p = 2
S = 5


###
int_const = mean_beta0[1]^2 - mean_beta0[1] * mean_beta0[2] + (1/3) * mean_beta0[2]^2 + 
  mean_beta1[1]^2 - mean_beta1[1] * mean_beta1[2] + (1/3) * mean_beta1[2]^2 - 
  2 * mean_beta0[1] * mean_beta1[1] - mean_beta0[1] * mean_beta1[2] - mean_beta0[2] * mean_beta1[1] - 
  (2/3) * mean_beta0[2] * mean_beta1[2]
int_const = sqrt(int_const)

# check, where mean_beta0[1] and mean_beta1[1] are both 0
#int_const2 = (1/3) * (mean_beta0[2] - mean_beta1[2])^2
#int_const2 = sqrt(int_const2)
# it worked

###
f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(Wass_dist)
}

unnormalized_f_dens_curve = function(x) f_dens(x,mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type)
integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = 0, upper = 1)$value


f_dens_curve = function(x) (1/int_const) * unnormalized_f_dens_curve(x)
f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)



###
hist(design_fast1$D[, S], probability = T, ylim = c(0, 2.5))
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)
hist(design_fast2$D[, S], probability = T, ylim = c(0, 2.5))
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)
hist(design_fast3$D[, S], probability = T, ylim = c(0, 2.5))
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)
hist(design_fast3.2$D[, S], probability = T, ylim = c(0, 2.5))
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)


