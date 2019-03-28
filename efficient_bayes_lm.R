############################################################
## Based on Fast Algorithm, Joseph et. al. 2018
############################################################

### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

## Like bayes_linear_regression.R, but I try to implement some of the techniques in
##  "Deterministic Sampling of Expensive Posteriors Using Minimum Energy Designs,"
##  Joseph 2018

library(transport)
library(mined)

source("smed_ms_functions.R")



### Testing

## Pick Parameters

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.05
var_e = 0.1 # same variance

N = 51 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5p = 103
K = 4 # ceiling(4* sqrt(p))
p = 1
xmin = 0
xmax = 1


## Running Algorithm


N = 11
source("smed_ms_functions.R")
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'efficient_bayes_FIXED_lm_N11.png')
#dev.off()


N = 51
source("smed_ms_functions.R")
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)

for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}

points(X_k, rep(0, N), col = 2)



N = 67
source("smed_ms_functions.R")
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)

for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}

points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_N67.png')
#dev.off()




# What if we draw from Uniform, instead of Lattice?
set.seed(1234)
N = 67
source("smed_ms_functions.R")
X_test = SMED_ms_fast2(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)

for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.007, i, col=4)
}

points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_unif_N67.png')
#dev.off()
