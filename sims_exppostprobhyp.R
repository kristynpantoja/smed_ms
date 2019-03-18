library(transport)
library(mined)

source("efficient_bayes_lm.R")

## --- Get the competing models and corresponding SMED design X --- ##

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.05
var_e = 0.1 # same variance

n = 51 # n = numCandidates should be largest prime number less than 100 + 5p = 103
K = 4 # ceiling(4* sqrt(p))
p = 1
xmin = 0
xmax = 1

#Warning message:
#In matrix(rep(D1, K), nrow = n, ncol = K) :
#  data length [188] is not a sub-multiple or multiple of the number of rows [51]


X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, xmin, xmax, K, p)

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

# just to see what it looks like again # nvm, too many points.
test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:n), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:n), col=4)
points(X_k, rep(0, n), col = 2)

## --- Simulate several Y and average to compare Posteriors H_i | X, Y --- ##

X = X_k

# suppose H_0 is true, i.e. y_i = x_i beta0 + epsilon_i
# then have:

numSims = 1000
Y = matrix(rep(NA, numSims * n), nrow = n, ncol = numSims) # each column is a new simulation of Y in R^n, for J = numSims total columns and n rows
for(i in 1:n){
  for(j in 1:numSims){
    Y[i, j] = rnorm(n = 1, mean = mean_beta0 * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))
  }
}

# just inspecting to see if they actually follow the line  y_i = x_i beta0
j = 40
plot(Y[ , j] ~ X)

# calculate expected evidence for each model, P(Y | H_i), i.e. averaged over beta (which is why we have mean_beta in the fmla)
# but it depends on x's in D... what do we do about those? some kind of a likelihood? Why would this work?




# calculate marginal y at each row and column, for each model
getMarginalY = function(x, mean_beta, i) dnorm(x, mean = mean_beta * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))

marginalY_H0 = matrix(rep(NA, numSims * n), nrow = n, ncol = numSims)
for(i in 1:n){
  for(j in 1:numSims){
    marginalY_H0[i, j] = getMarginalY(Y[i, j], mean_beta0, i)
  }
}

marginalY_H1 = matrix(rep(NA, numSims * n), nrow = n, ncol = numSims)
for(i in 1:n){
  for(j in 1:numSims){
    marginalY_H1[i, j] = getMarginalY(Y[i, j], mean_beta1, i)
  }
}
# calculate expected marginal y for each model by averaging
expected_marginalY_H0 = mean(marginalY_H0)
expected_marginalY_H1 = mean(marginalY_H1)

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (0.5 * expected_marginalY_H0 + 0.5 * expected_marginalY_H1)
expected_post_H1 = expected_marginalY_H1 / (0.5 * expected_marginalY_H0 + 0.5 * expected_marginalY_H1)
