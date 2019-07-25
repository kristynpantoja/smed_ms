###### playing around with scatterplot3d ######

library(scatterplot3d)
## example 1
z <- c(2, 2)
x <- c(1, 1)
y <- c(0, 1)
scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), 
              highlight.3d=TRUE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d")
## example 2
N = 10
z <- rep(2, N * N)
x <- rep(seq(from = 0, to = 1, length.out = N), each = N)
y <- rep(seq(from = 0, to = 1, length.out = N), times = N)
scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), 
              highlight.3d=TRUE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d")

## example 3
N = c(10, 12)
z <- rep(2, N[1] * N[2])
x <- rep(seq(from = 0, to = 1, length.out = N[1]), each = N[2])
y <- rep(seq(from = 0, to = 1, length.out = N[2]), times = N[1])
scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), 
              highlight.3d=TRUE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d")

## see candidates
z <- rep(0.5, 100 * 100)
x <- candidates[ , 1]
y <- candidates[ , 2]
scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), 
              highlight.3d=TRUE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d")




###### some 2D fns to compare ######

## example 1
f0 = function(x, y) return(0)
f1 = function(x, y) return(0.5)


## example 2
f0 = function(x, y) return(0)
f1 = function(x, y) return(1 - 0.5 * (x + y))



#####################
###### 2-DIM'L ######
#####################

# Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly_2d = function(x, var_mean, var_e, type, var_margy){
  # type:
  #   1 for linear model without slope
  #   2 for linear model with slope
  #   3 for quadratic model with slope
  if(!is.null(type)){
    if(type == 4) var_mean[1] + x[1]^2 * var_mean[2] + x[2]^2 * var_mean[3] + var_e
    else stop(paste("invalid type given : ", type))
  } else{
    if(!is.null(var_margy)) var_margy(x = x, var_mean = var_mean, var_e = var_e)
    else stop("cannot compute var_marginaly: no type given, and no var_margy fn given")
  }
}

# charge function evaluated at x
q_2d = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly_2d(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly_2d(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(1.0 / Wass_dist^(1 / (2 * p)))
}


###################
### one at time ###
###################

f_min_2d = function(candidate, D, k, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                 f0, f1, type, var_margy0, var_margy1, p, log_space = FALSE){
  if(log_space == FALSE) {
    
    result = q_2d(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)^k * sum(apply(D, 1, function(x_i) (q_2d(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p) / sqrt(sum((x_i - candidate)^2)))^k))
    return(result)
  } else{
    # if has logSumExp library
    terms = sapply(D, function(x_i) k * log(q(candidate, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) + k * log(q(x_i, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0, var_margy1, p)) - k * log(sqrt((x_i - candidate)^2)))
    result = exp(logSumExp(terms))
    return(result)
  }
}

MED_ms_2d = function(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                  f0 = NULL, f1 = NULL, type = NULL, var_margy0 = NULL, var_margy1 = NULL, 
                  N = 11, numCandidates = 10^5, k = 4, p = 2, xmin = 0, xmax = 1, log_space = FALSE, 
                  genCandidates = 1, initialpt = 1){
  # var_margy0 and var_margy1 : functions that take in x, var_mean, var_e
  
  if(is.null(type) & is.null(f0) & is.null(f1) & is.null(var_margy0) & is.null(var_margy0)) stop("must specify model type and/or model")
  
  # Create hypothesized models
  if(is.null(f0)){
    if(type[1] == 4) f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2]
    else stop("type[1] is invalid and f0 is not provided")
  }
  if(is.null(f1)){
    if(type[1] == 4) f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2]
    else stop("type[2] is invalid and f1 is not provided")
  }
  
  # -- Generate Candidate Points -- #
  if(length(numCandidates) == 1) numCandidates = c(floor(sqrt(numCandidates)), floor(sqrt(numCandidates)))
  #### CHECK ###
  if(genCandidates == 1){
    candidates_x1 = seq(from = xmin, to = xmax, length.out = numCandidates[1])
    candidates_x2 = seq(from = xmin, to = xmax, length.out = numCandidates[2])
    candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
                       rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
  }
  if(genCandidates == 2){
    candidates_x1 = sort(runif(numCandidates[1], min = xmin, max = xmax))
    candidates_x2 = sort(runif(numCandidates[2], min = xmin, max = xmax))
    candidates = cbind(rep(candidates_x1, each = numCandidates[2]), 
                       rep(candidates_x2, times = numCandidates[1])) # each row is a candidate (x1, x2)
  }
  
  # -- Initialize 1st Design Point in D -- #
  D = matrix(rep(NA, N * 2), N, 2)
  if(initialpt == 2){
    f0_cand = apply(candidates, 1, FUN = f0)
    f1_cand = apply(candidates, 1, FUN = f1)
    xinitind = which.max(abs(f0_cand - f1_cand))
    D[1, ] = candidates[xinitind, ]
  } else{
    D[1, ] = optimize(function(x) q_2d(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, 
                                  type, var_margy0, var_margy1, p), interval = c(xmin, xmax))$minimum
  }
  
  for(i in 2:N){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_2d(x, D[1:(i - 1), , drop = FALSE], k, mean_beta0, mean_beta1, 
                                                            var_mean0, var_mean1, var_e, f0, f1, 
                                                            type, var_margy0, var_margy1, p))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt, ]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
  }
  return(D)
}




######## test it out!

library(expm)
library(matrixStats)
source("/Users/kristyn/Documents/research/smed_ms/med_ms_functions.R")


# Model 1 : f0(x, y) = 0, f1(x, y) = 0.5
mean_beta0 = c(0, 0, 0) # null model
mean_beta1 = c(0.5, 0, 0) # alternative model
var_mean0 = diag(c(0.005, 0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2] # null
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2] # alt
N = 25
type = c(4, 4)
p = 3
# for one-at-a-time algorithm:
numCandidates = 10e3 # suggested 10^5
k = 4
#design = MED_ms_2d(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                     var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                     f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                     N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
#                     genCandidates = 1, initialpt = 1)
## see design points
z <- rep(0.5, N)
x <- design[ , 1]
y <- design[ , 2]
scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), 
              highlight.3d=TRUE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d")
#hists
hist(design[ , 1])
hist(design[ , 2])



# Model 2 : f0(x, y) = 0, f1(x, y) = 1 - 0.5 * (x + y)
mean_beta0 = c(0, 0, 0) # null model
mean_beta1 = c(1, -0.5, -0.5) # alternative model
var_mean0 = diag(c(0.005, 0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2] # null
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2] # alt
N = 25
type = c(4, 4)
p = 3
# for one-at-a-time algorithm:
numCandidates = 10e3 # suggested 10^5
k = 4
#design2 = MED_ms_2d(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                   var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                   f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                   N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
#                   genCandidates = 1, initialpt = 1)
## see design points
# function 1 - 0.5 * (x + y)
numPts = 25
xy <- cbind(rep(seq(from = 0, to = 1, length.out = numPts), each = numPts), 
            rep(seq(from = 0, to = 1, length.out = numPts), times = numPts))
f1_eval_xy = apply(xy, 1, function(x) 1 - 0.5 * x[1] - 0.5 * x[2])
z <- f1_eval_xy
x <- xy[ , 1]
y <- xy[ , 2]
sc = scatterplot3d(x, y, z, 
              xlim = c(0, 1), ylim = c(0, 1), color = 2,
              highlight.3d = FALSE, col.axis="blue", grid = TRUE,
              col.grid="lightblue", main="scatterplot3d", pch = ".")
# points
f1_eval_design2 = apply(design2, 1, function(x) 1 - 0.5 * x[1] - 0.5 * x[2])
z <- f1_eval_design2
x <- design2[ , 1]
y <- design2[ , 2]
sc$points3d(x, y, z, xlim = c(0, 1), ylim = c(0, 1))
# hists
hist(design2[ , 1])
hist(design2[ , 2])


# Model 2.1 : f0(x, y) = 0, f1(x, y) = 1 - 0.5 * (x + y) MORE POINTS
mean_beta0 = c(0, 0, 0) # null model
mean_beta1 = c(1, -0.5, -0.5) # alternative model
var_mean0 = diag(c(0.005, 0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2] # null
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[2] # alt
N = 50
type = c(4, 4)
p = 3
# for one-at-a-time algorithm:
numCandidates = 10e3 # suggested 10^5
k = 4
#design2.1 = MED_ms_2d(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#                    var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
#                    f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
#                    N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
#                    genCandidates = 1, initialpt = 1)
## see design points
# function 1 - 0.5 * (x + y)
numPts = 25
xy <- cbind(rep(seq(from = 0, to = 1, length.out = numPts), each = numPts), 
            rep(seq(from = 0, to = 1, length.out = numPts), times = numPts))
f1_eval_xy = apply(xy, 1, function(x) 1 - 0.5 * x[1] - 0.5 * x[2])
z <- f1_eval_xy
x <- xy[ , 1]
y <- xy[ , 2]
sc = scatterplot3d(x, y, z, 
                   xlim = c(0, 1), ylim = c(0, 1), color = 2,
                   highlight.3d = FALSE, col.axis="blue", grid = TRUE,
                   col.grid="lightblue", main="scatterplot3d", pch = ".")
# points
f1_eval_design2.1 = apply(design2.1, 1, function(x) 1 - 0.5 * x[1] - 0.5 * x[2])
z <- f1_eval_design2.1
x <- design2.1[ , 1]
y <- design2.1[ , 2]
sc$points3d(x, y, z, xlim = c(0, 1), ylim = c(0, 1))
# hists
hist(design2.1[ , 1])
hist(design2.1[ , 2])



