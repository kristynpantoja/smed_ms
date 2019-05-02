##########################################
### --- evaluating both algorithms --- ###
##########################################

setwd("/Users/kristyn/Documents/research/smed_ms")
library(transport)
library(mined)
library(expm)
source("smed_ms_functions.R")

##########################
### --- parameters --- ###
##########################

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 0.01 # same variance

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

N = 67

# for fast algorithm:
K = 40 # ceiling(4* sqrt(p))
p = 1 * 2

# for one-at-a-time algorithm:
numCandidates = 10^5 # suggested 10^5
k = 4



###################################
### --- efficient algorithm --- ###
###################################

# TIME IT!
set.seed(1)
ptm <- proc.time()
X_test_efficient = SMED_ms_fast(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                                var_e = var_e, var_mean = var_mean, N = N, 
                                xmin = xmin, xmax = xmax, K = K, p = 2)
proc.time() - ptm # elapsed = 581.156

# HISTOGRAM.
test_k = K
hist(X_test_efficient$D[ ,test_k], freq = T)
mean(X_test_efficient$D[ ,test_k]) # 0.6143494 (compare to 0.6129598 in K = 100)
sd(X_test_efficient$D[ ,test_k]) # 0.2557405 (compare to 0.2540533 in K = 100)
#dev.copy(png, 'fast_N67_K40_k40_p2_hist.png')
#dev.off()


# PICTURES.
test_k = K
#X_k = sort(X_test_efficient$D[ ,test_k])
#curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
#curve(f1, col = 1, add = TRUE)
#text(X_test_efficient$D[ ,test_k], f0(X_test_efficient$D[ ,test_k]), c(1:N), col=4)
#text(X_test_efficient$D[ ,test_k], f1(X_test_efficient$D[ ,test_k]), c(1:N), col=4)
#points(X_k, rep(0, 11), col = 2)


# see what the last design looks like
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

#dev.copy(png, "fast_pattern2_N67_K40_k40_p2.png")
#dev.off()


### --- Computing Posterior Probabilities of Hypotheses for N = 67, K = 40, p = 2 --- ###
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
expected_marginalY_H0 = mean(marginalY_H0) # 2.756252
expected_marginalY_H1 = mean(marginalY_H1) # 0.6076236

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.819368
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.180632

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 4.536117 > 1, Null hypothesis is supported


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
expected_marginalY_H0 = mean(marginalY_H0) # 0.617338
expected_marginalY_H1 = mean(marginalY_H1) # 2.764161

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.1825634
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.8174366

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 0.2233365 < 1, Null hypothesis is not supported


### --- Distances between x's and y's N = 67 --- ###

Xsort = sort(X_test_efficient$D[ ,test_k])
Y0 = f0(Xsort)
Y1 = f1(Xsort)
diffX = Xsort[-1] - Xsort[-length(Xsort)]
diffY = Y0 - Y1

sum(diffX) # 0.9599299
sum(diffY) # 20.5807

hist(diffX) # right-skewed, like dist of X
hist(diffY) # left-skewed, but way less severe drop off

mean(diffX) # 0.01454439
mean(diffY) # 0.3071747

sd(diffX) # 0.007635743
sd(diffY) # 0.1278703



### --- Standardized Distances between x's and y's ###


### --- Compute Criterion ###

totalPE = function(D, N, mean_beta0, mean_beta1, var_e, var_mean){
  if(N != length(D)) stop("N is not the same as length of D")
  numPairs = N * (N - 1) / 2
  pairwise_PEs = rep(NA, numPairs)
  counter = 1
  qD = sapply(FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean), D)
  for(i in 1:(N - 1)){
    for(j in (i + 1):N){
      pairwise_PEs[counter] = qD[i] * qD[j] / (D[i] - D[j])^2
      counter = counter + 1
    }
  }
  return(sum(pairwise_PEs))
}

totalPE_1atatime = totalPE(X_test_1atatime, N, mean_beta0, mean_beta1, var_e, var_mean)
totalPE_each_k = rep(NA, K)
for(k in 1:K){
  totalPE_each_k[k] = totalPE(X_test_efficient$D[ , k], N, mean_beta0, mean_beta1, var_e, var_mean)
}
plot(totalPE_each_k, type = "l")
abline(a = totalPE_1atatime, b = 0, col = 2)

?plot
#######################################
### --- one-at-a-time algorithm --- ###
#######################################


## Running Algorithm

# TIME IT!
set.seed(1)
ptm <- proc.time()
X_test_1atatime = SMED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1,
                          var_e = var_e, var_mean = var_mean, N = N, numCandidates = numCandidates, 
                          k = k, xmin = xmin, xmax = xmax)
proc.time() - ptm # 1022.885

# HISTOGRAM.
hist(X_test_1atatime, freq = T)
mean(X_test_1atatime) # 0.6882634
sd(X_test_1atatime) # 0.2182556

#dev.copy(png, 'oneatatime_seq_pattern2_N67_nC10to5_hist.png')
#dev.off()

### --- Distances between x's and y's N = 67 --- ###

Xsort = sort(X_test_1atatime)
Y0 = f0(Xsort)
Y1 = f1(Xsort)
diffX = Xsort[-1] - Xsort[-length(Xsort)]
diffY = Y0 - Y1

sum(diffX) # 0.8859489
sum(diffY) # 23.05683

hist(diffX) # right-skewed, like dist of X
hist(diffY) # left-skewed, but less severe drop off

mean(diffX) # 0.01342347
mean(diffY) # 0.3441317

sd(diffX) # 0.01187598
sd(diffY) # 0.1091278


# PICTURES.

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 3))
curve(f1, col = 1, add = TRUE)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_test_1atatime[i], y[i], i, col=4)
}
points(X_test_1atatime, rep(0, n), col = 2, pch = "*")
lines(X_test_1atatime, y, col = 3)

#dev.copy(png, 'oneatatime_seq_pattern2_N67_nC10to5.png')
#dev.off()


### --- Computing Posterior Probabilities of Hypotheses for One-at-a- Time --- ###
# Simulate several Y and average to compare Posteriors H_i | X, Y

X = X_test_1atatime$D[ ,test_k]

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
expected_marginalY_H0 = mean(marginalY_H0) # 2.756252
expected_marginalY_H1 = mean(marginalY_H1) # 0.6076236

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.819368
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.180632

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 4.536117 > 1, Null hypothesis is supported


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
expected_marginalY_H0 = mean(marginalY_H0) # 0.617338
expected_marginalY_H1 = mean(marginalY_H1) # 2.764161

# calculate expected posterior probability of each model (equal prior on models)
expected_post_H0 = expected_marginalY_H0 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.1825634
expected_post_H1 = expected_marginalY_H1 / (expected_marginalY_H0 + expected_marginalY_H1) # 0.8174366

BayesFactor_01 = expected_marginalY_H0 / expected_marginalY_H1 # since 0.2233365 < 1, Null hypothesis is not supported







