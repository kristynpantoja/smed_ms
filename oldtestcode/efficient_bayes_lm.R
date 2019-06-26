#################################################
## Based on Fast Algorithm, Joseph et. al. 2018 #
#################################################

### Linear Regression Different Slopes (both intercept at 0, same error variance) 

## Like bayes_linear_regression.R, but I try to implement some of the techniques in
##  "Deterministic Sampling of Expensive Posteriors Using Minimum Energy Designs,"
##  Joseph 2018



### --- Libraries / Sourced Scripts --- ###
library(transport)
library(mined)
library(expm)

source("smed_ms_functions.R")













###############
### Testing ###
###############


### --- Pick Parameters --- ###

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 0.01 # same variance

p = 1 
xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

# parameters we change sometimes: 
N = 51 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5 * p = 103
K = 4 # ceiling(4* sqrt(p))













############################
## Running Algorithm, K = 4
############################
K = 4

### --- N = 11 --- ###

N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_N11.png')
#dev.off()



### --- N = 67 --- ###
# does not vary by much

N = 67
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N67.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_N67.png')
#dev.off()


### --- N = 67, but with uniformly selected candidates --- ###


# What if we use uniform, instead of Lattice?
# high variability!!! why??? because numCandidates at each step is equal to N

N = 67
X_test = SMED_ms_fast2(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastunif_N67.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastunif_pattern_N67.png')
#dev.off()













#############################
## Running Algorithm, K = 20
#############################
K = 20

### --- N = 11 --- ###

N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11_K20.png')
#dev.off()



### --- N = 67 --- ###

N = 67

# TIME IT!
ptm <- proc.time()
X_test_K20 = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
proc.time() - ptm # elapsed = 3507.396

test_k = 20

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 3))
curve(f1, col = 1, add = TRUE)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_test_K20$D[i ,test_k], y[i], i, col=4)
}
points(X_test_K20$D[ ,test_k], rep(0, N), col = 2, pch = "*")
lines(X_test_K20$D[ ,test_k], y, col = 3)

#dev.copy(png, 'fast_pattern2_N67_K20_k20.png')
#dev.off()













#############################
## Running Algorithm, K = 100
#############################
K = 100

### --- N = 11 --- ###

N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11_K100.png')
#dev.off()



### --- N = 67 --- ###

N = 67

# TIME IT!
ptm <- proc.time()
X_test_efficient = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
proc.time() - ptm # elapsed = 3507.396

# HISTOGRAM.
test_k = 100
hist(X_test_efficient$D[ ,test_k], freq = T)
mean(X_test_efficient$D[ ,test_k]) # 0.6129598
sd(X_test_efficient$D[ ,test_k]) # 0.2540533
#dev.copy(png, 'fast_pattern2_N67_K100_k100_2_hist.png')
#dev.off()


# PICTURES.
#test_k = 100
#X_k = sort(X_test_efficient$D[ ,test_k])
#curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
#curve(f1, col = 1, add = TRUE)
#text(X_test_efficient$D[ ,test_k], f0(X_test_efficient$D[ ,test_k]), c(1:N), col=4)
#text(X_test_efficient$D[ ,test_k], f1(X_test_efficient$D[ ,test_k]), c(1:N), col=4)
#points(X_k, rep(0, N), col = 2)


# see what the last design looks like
test_k = 100
X_k = sort(X_test_efficient$D[ ,test_k])
#curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0.0, 2.5))
#curve(f1, col = 1, add = TRUE)
#for(i in 1:N){
#  text(X_test_efficient$D[i ,test_k], f1(X_test_efficient$D[i ,test_k]) + i * 0.03, i, col=4)
#}
#points(X_k, rep(0, N), col = 2)

#ks = seq(from = 0, to = 100, by = 10)
#ks[1] = 1
#for(i in 1:length(ks)){
#  name = paste("fast_pattern_N67_K100_k", ks[i], ".png", sep = "")
#  
#  test_k = ks[i]
#  X_k = sort(X_test_efficient$D[ ,test_k])
#  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", 
#        ylim = c(0.0, 2.5))
#  curve(f1, col = 1, add = TRUE)
#  for(i in 1:N){
#    text(X_test_efficient$D[i ,test_k], f1(X_test_efficient$D[i ,test_k]) + i * 0.03, i, col=4)
#  }
#  points(X_k, rep(0, N), col = 2)
#  dev.copy(png, name)
#  dev.off()
#}

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 3))
curve(f1, col = 1, add = TRUE)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_test_efficient$D[i ,test_k], y[i], i, col=4)
}
points(X_test_efficient$D[ ,test_k], rep(0, N), col = 2, pch = "*")
lines(X_test_efficient$D[ ,test_k], y, col = 3)

#dev.copy(png, "fast_pattern2_N67_K100_k100_2.png")
#dev.off()


### --- Computing Posterior Probabilities of Hypotheses for N = 67 --- ###
# Simulate several Y and average to compare Posteriors H_i | X, Y

X = X_test_efficient$D[ ,test_k]

## 1st suppose H_0 is true, i.e. y_i = x_i beta0 + epsilon_i
# then have:

numSims = 1000
# each column is a new simulation of Y in R^n, for J = numSims total columns and n rows
Y = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims)
for(i in 1:N){
  for(j in 1:numSims){
    Y[i, j] = rnorm(n = 1, mean = mean_beta0 * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))
  }
}

# just inspecting to see if they actually follow the line  y_i = x_i beta0
randSim = sample(1:numSims, 1); plot(Y[ , randSim] ~ X)

# calculate expected evidence for each model, P(Y | H_i), 
# i.e. averaged over beta (which is why we have mean_beta in the fmla)

# calculate marginal y (given H0) at each row and column, for each model
getMarginalY = function(x, mean_beta, i) dnorm(x, mean = mean_beta * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))

marginalY_H0 = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims)
for(i in 1:N){
  for(j in 1:numSims){
    marginalY_H0[i, j] = getMarginalY(Y[i, j], mean_beta0, i)
  }
}

marginalY_H1 = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims)
for(i in 1:N){
  for(j in 1:numSims){
    marginalY_H1[i, j] = getMarginalY(Y[i, j], mean_beta1, i)
  }
}

# calculate expected marginal y for each model by averaging
expected_marginalY_H0 = mean(marginalY_H0) # 2.759501
expected_marginalY_H1 = mean(marginalY_H1) # 0.6073309

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.8196135
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.1803865

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 4.543653 > 1, Null hypothesis is supported


## 2nd suppose H_1 is true, i.e. y_i = x_i beta1 + epsilon_i
# then have:

numSims = 1000
Y = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims) # each column is a new simulation of Y in R^n, for J = numSims total columns and n rows
for(i in 1:N){
  for(j in 1:numSims){
    Y[i, j] = rnorm(n = 1, mean = mean_beta1 * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))
  }
}

# just inspecting to see if they actually follow the line  y_i = x_i beta1
randSim = sample(1:numSims, 1); plot(Y[ , randSim] ~ X)

# calculate expected evidence for each model, P(Y | H_i), 
# i.e. averaged over beta (which is why we have mean_beta in the fmla)

# calculate marginal y at each row and column, for each model
getMarginalY = function(x, mean_beta, i) dnorm(x, mean = mean_beta * X[i], sd = sqrt(var_e + X[i]^2 * var_mean))

marginalY_H0 = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims)
for(i in 1:N){
  for(j in 1:numSims){
    marginalY_H0[i, j] = getMarginalY(Y[i, j], mean_beta0, i)
  }
}

marginalY_H1 = matrix(rep(NA, numSims * N), nrow = N, ncol = numSims)
for(i in 1:N){
  for(j in 1:numSims){
    marginalY_H1[i, j] = getMarginalY(Y[i, j], mean_beta1, i)
  }
}
# calculate expected marginal y for each model by averaging
expected_marginalY_H0 = mean(marginalY_H0) # 0.6139273
expected_marginalY_H1 = mean(marginalY_H1) # 2.761371

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.1818883
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.8181117

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 0.222327 < 1, Null hypothesis is not supported


### --- Distances between x's and y's N = 67 --- ###

Xsort = sort(X_test_efficient$D[ ,test_k])
Y0 = f0(Xsort)
Y1 = f1(Xsort)
diffX = Xsort[-1] - Xsort[-length(Xsort)]
diffY = Y0 - Y1

sum(diffX) # 0.9599863
sum(diffY) # 20.53415

hist(diffX) # right-skewed, like dist of X
hist(diffY) # left-skewed, but less severe drop off

mean(diffX) # 0.01454525
mean(diffY) # 0.3064799

sd(diffX) # 0.00766221
sd(diffY) # 0.1270267


















#############################################
## Running Algorithm, K = 40, Lattice Version
#############################################
# just because that's the function they use in their code
K = 100

### --- N = 11 --- ###

N = 11
X_test = SMED_ms_fast3(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastlattice_N11_K100.png')
#dev.off()


### --- N = 67 --- ###

N = 67
X_test = SMED_ms_fast3(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = K

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0.0, 2.5))
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
}
points(X_k, rep(0, N), col = 2)

# looks similar to seq version

ks = seq(from = 0, to = 100, by = 10)
ks[1] = 1
for(i in 1:length(ks)){
  name = paste("fast_pattern_N67_K100_k", ks[i], ".png", sep = "")
  
  test_k = ks[i]
  X_k = sort(X_test$D[ ,test_k])
  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", 
        ylim = c(0.0, 2.5))
  curve(f1, col = 1, add = TRUE)
  for(i in 1:N){
    text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
  }
  points(X_k, rep(0, N), col = 2)
  dev.copy(png, name)
  dev.off()
}



















