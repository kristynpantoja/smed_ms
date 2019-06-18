library(transport)
library(mined)

source("efficient_bayes_lm.R")

## --- Get the competing models and corresponding SMED design X --- ##

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 0.01 # same variance

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

#dev.copy(png,'sim_H0')
#dev.off()

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
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1)
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1)

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since > 1, Null hypothesis is supported
















# suppose H_1 is true, i.e. y_i = x_i beta0 + epsilon_i
# then have:

numSims = 1000
Y = matrix(rep(NA, numSims * n), nrow = n, ncol = numSims) # each column is a new simulation of Y in R^n, for J = numSims total columns and n rows
for(i in 1:n){
  for(j in 1:numSims){
    Y[i, j] = rnorm(n = 1, mean = mean_beta1 * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))
  }
}

# just inspecting to see if they actually follow the line  y_i = x_i beta0
j = 40
plot(Y[ , j] ~ X)

#dev.copy(png,'sim_H1')
#dev.off()


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
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1)
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1)

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since > 1, Null hypothesis is supported






















########################################################################################
################################### TRY AGAIN ##########################################
########################################################################################
mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.1; var_mean1 = var_mean0; # variance on beta
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

N = 67
type = c(1, 1)
p = 1

# for one-at-a-time algorithm:
numCandidates = 10^3 # suggested 10^5
k = 4
# One-at-a-Time Algorithm
X_1atatime = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                    var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                    f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                    N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
                    genCandidates = 1, initialpt = 1)


library(mvtnorm)
numSims = 1000

model_evidence = function(Y, X, mean_beta, var_mean, var_e){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  N = length(Y)
  marginaly_mean = X %*% mean_beta
  marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

simulateY = function(X, mean_beta, var_mean, var_e, numSims, plotrandsim = FALSE){
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    beta = rnorm(n = 1, mean = mean_beta, sd = sqrt(var_mean))
    for(i in 1:N){
      Y[i, j] = rnorm(n = 1, mean = beta * X[i], sd = sqrt(var_e))
    }
  }
  if(plotrandsim == TRUE){
    randSim = sample(1:numSims, 1)
    plot(Y[ , randSim] ~ X)
  }
  return(Y)
}

# testing - it works fine! :D

X = as.matrix(X_1)
X = as.matrix(X_1atatime)

simY = simulateY(X, mean_beta = mean_beta0, var_mean0, var_e, numSims)
randSim = sample(1:numSims, 1); plot(simY[ , 104] ~ X, xlim = c(0, 1), ylim = c(0, 1))

simPostH0 = rep(NA, numSims)
simPostH1 = rep(NA, numSims)
simBF01 = rep(NA, numSims)

j = 104

for(j in 1:numSims){
  Y = simY[, j]
  # get model evidences
  simEvidenceH0 = model_evidence(Y, X, mean_beta0, var_mean0, var_e)
  simEvidenceH1 = model_evidence(Y, X, mean_beta1, var_mean1, var_e)
  # calculate posterior probabilities of models
  simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
  simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
  # calculate bayes factor
  simBF01[j] = simPostH0[j] / simPostH1[j]
}
expected_postH0 = mean(simPostH0)
expected_postH1 = mean(simPostH1)
expected_BF01 = mean(simBF01)


