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

n = 11 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5p = 103
K = 4 # ceiling(4* sqrt(p))
p = 1
xmin = 0
xmax = 1

## Running Algorithm

X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, n, xmin, xmax, K, p)

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_k, rep(0, n), col = 2)



## Add errors for marginal of y | H_i
X_test_errors = sapply(X_k, function(x) var_marginaly(x, var_e, var_mean))
polygon(x = c(X_k, rev(X_k)),y = c(f0(X_k) - 2 * X_test_errors, rev(f0(X_k) + 2 * X_test_errors)),
        col =  adjustcolor("cyan", alpha.f = 0.10), border = NA)
polygon(x = c(X_k, rev(X_k)),y = c(f1(X_k) - 2 * X_test_errors, rev(f1(X_k) + 2 * X_test_errors)),
        col =  adjustcolor("pink", alpha.f = 0.10), border = NA)



## Debugging to find out candidates calculations
test_k = 1
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_test$D[ ,test_k], rep(0, n), col=2)

length(X_test$candidates[[3]]) # why is it 132?


