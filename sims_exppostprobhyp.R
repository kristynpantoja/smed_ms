
mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.005; var_mean1 = var_mean0; # variance on beta
var_e = 0.025 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

N = 67
type = c(1, 1)
p = 1


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
D_space = seq(xmin, xmax, length.out = N)
X = as.matrix(D_space)

simY = simulateY(X, mean_beta = mean_beta0, var_mean0, var_e, numSims)
randSim = sample(1:numSims, 1); plot(simY[ , randSim] ~ X, xlim = c(0, 1), ylim = c(0, 1))

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


