# testing higher dimensions

# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

# --- libraries --- #
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# --- sources to generate MEDs --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))

source(paste(functions_home, "/simulate_seqMED_multidim.R", sep = ""))

##########
# CASE 1 #
##########

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

# MED design #
typeT = 3
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
xmin = -1
xmax = 1

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
curve(f0, col = 2, lwd = 5, xlim = c(xmin, xmax), xlab = "x", ylab = "f")
curve(f1, add = T, col = 5, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)
legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))

# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^3
k = 4
xmin = -1
xmax = 1
p = 1
N = 15

numSeq = 3
N_seq = 5
alpha_seq = 1
smmed1 = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
                         numCandidates, k, xmin, xmax, p, numSeq = numSeq, 
                         N_seq = N_seq, alpha_seq = alpha_seq,seed = 12)
plot(x = 1:N, y = smmed1$D, type = "l")



smmed2 = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
                         numCandidates, k, xmin, xmax, p, numSeq = numSeq, 
                         N_seq = N_seq, alpha_seq = alpha_seq,seed = 12)
plot(x = 1:N, y = smmed2$D, type = "l")

################################
## attempt at multidimensional #
################################

# simulate_seqMED_multidim args
true_beta = betaT
true_type = typeT
mean_beta0 = mu0
mean_beta1 = mu1
var_beta0 = V0
var_beta1 = V1
var_e = sigmasq
# f0 = NULL
# f1 = NULL
type = type01
# numCandidates = 10^5
# k = 4, 
# xmin = 0
# xmax = 1
# p = 1
initD = rnorm(5)
inity = f1(initD) + rnorm(5, 0, sqrt(sigmasq))
numSeq = 3 - 1
N_seq = 5
alpha_seq = 1
buffer_seq = 0
wasserstein0 = 1
genCandidates = 1
seed = 12


# # add_MED_ms_oneatatime_data_multidim args
# initD = D # D gets updated
# y = y # y also gets updated
# N2 = N_seq[t]
# numCandidates = numCandidates
# alpha = alpha_seq[t]
# buffer = buffer_seq[t]














########################################
###########################################
##############################################

# testing higher dimensions

# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

# --- libraries --- #
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# --- sources to generate MEDs --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))

source(paste(functions_home, "/simulate_seqMED_multidim.R", sep = ""))

##########
# CASE 1 #
##########
mu_full = c(0, 0, 0)


indices0 = c(1, 2)
indices1 = 1:length(mu_full)
mu0 = mu_full[indices0]
mu1 = mu_full[indicse1]

betaT = c(10, 11, 12)
indicesT = indices0
sigmasq01 = 0.25
sigmasq = 0.1

# MED design #
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]
xmin = -1
xmax = 1

fT = function(x) betaT %*% x[indicesT]

# checks
indicesT <= max(indices0, indices1)
length(betaT) == length(indicesT)

# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
numCandidates = 10^3
k = 4
xmin = -1
xmax = 1
p = 1
N = 15

numSeq = 3
N_seq = 5
alpha_seq = 1


################################
## multidim settings in slides #
################################

# libraries and helper fns
library(mvtnorm)
invisible(library(glmnet))
standardizeXY = function(X, Y){
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Center Y
  Ymean = mean(Y)
  Ytilde = Y - mean(Y)
  
  # Center and scale X
  Xmeans = colMeans(X)
  Xcen = X - matrix(Xmeans, n, p, byrow=T)
  normsX2 = colSums(Xcen^2) / n
  weights = 1 / sqrt(normsX2) # should weights be the vector, or the matix?
  Xtilde = Xcen %*% diag(weights)
  
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# what is the saturated model? has all coefficients, even the useless ones
# what is the sparse model? has only the useful coefficients (all coefficients except the useless ones)

# generate beta_sat
set.seed(12)
numUseful = 5 # number of useful parameters (the ones used to generate the data)
numUseless = 10 # number of useless parameters
numSat = numUseful + numUseless

# get the true parameters (and their indices)
beta_useful = rep(1, numUseful) * sample(c(-1, 1), numUseful, replace = TRUE) #runif(n = numUseful, min = 1.25, max = 1.5) * sample(c(-1, 1), numUseful, replace = TRUE)
beta_useless = rep(0, numUseless)
indices = sample(1:numSat, numSat, replace = FALSE)
ind_useful = indices[1:numUseful]
ind_useless = indices[(numUseful + 1):numSat]
beta_sat = rep(NA, numSat)
beta_sat[ind_useful] = beta_useful
beta_sat[-ind_useful] = beta_useless
# beta_sat = sample(c(beta_useful, beta_useless), p_sat, replace = FALSE)
plot(x = 1:numSat, y = beta_sat); abline(h = 0, lty = 3)

# define true model
beta0 = 0
truef = function(x) beta0 + t(beta_useful) %*% x

# generate data
set.seed(123*sample(1:50000, 1))
xmin = -1
xmax = 1
n_obs = 50 # sample size
sigma = 1 # noise standard deviation

Sigma_sat = matrix(0.5, numSat, numSat) + diag(rep(1-0.5, numSat)) # covariance for design X
X_sat = rmvnorm(n_obs, mean = rep(0, numSat), sigma = Sigma_sat) # matrix of covariates
X_sat = as.data.frame(X_sat)

Sigma_useful = Sigma_sat[ind_useful, ind_useful]
X_useful = X_sat[ , ind_useful]

Y = beta0 + as.matrix(X_useful) %*% beta_useful + sigma * rnorm(n_obs) # response from linear model

# creat train and test sets
proportion_split = 0.5
train_index = sample(1:n_obs, floor(n_obs * proportion_split))
Xtrain = X_sat[train_index, ]
Xtest = X_sat[-train_index, ]
Ytrain = Y[train_index]
Ytest = Y[-train_index]

trainset = data.frame(Xtrain, "response" = Ytrain)
testset = data.frame(Xtest, "response" = Ytest)
response_index = dim(trainset)[2]

# train and test set for cv?

# fit saturated linear model
lmsat = lm(response ~ . + 0, data = trainset)
Yhat_lmsat = predict(lmsat, newdata = trainset)
# all.equal(as.vector(Yhat_lmsat), as.vector(as.matrix(Xtrain) %*% lmsat$coefficients))
rsq_lmsat = cor(Ytrain, Yhat_lmsat)^2
rss_lmsat = sum((Ytrain - Yhat_lmsat)^2)
Yhat_test_lmsat = predict(lmsat, newdata = testset)
rss_test_lmsat = sum((Ytest - Yhat_test_lmsat)^2)

# fit sparse model (no useless variables)
lmsparse = lm(response ~ . + 0, data = trainset[ , c(ind_useful, response_index)])
Yhat_lmsparse = predict(lmsparse, newdata = trainset) # note: same as if put trainset[ , ind_useful], because it just grabs the variables that can be used in the model
# all.equal(as.vector(Yhat_lmsparse), as.vector(as.matrix(Xtrain[ , ind_useful]) %*% lmsparse$coefficients))
rsq_lmsparse = cor(Ytrain, Yhat_lmsparse)^2
rss_lmsparse = sum((Ytrain - Yhat_lmsparse)^2)
Yhat_test_lmsparse = predict(lmsparse, newdata = testset)
rss_test_lmsparse = sum((Ytest - Yhat_test_lmsparse)^2)

# fit ridge regression
# get lambda_seq
n_train = dim(Xtrain)[1]
p = dim(Xtrain)[2]
# XYstd = standardizeXY(X_sat, Y)
# Xtilde = XYstd$Xtilde
# Ytilde = XYstd$Ytilde
n_lambda = 30
lambda_max = 1 # for lasso: max(abs(crossprod(Xtilde, Ytilde) / n))
lambda_seq <- 10^seq(-3, 5, length.out = 100)
# lambda_seq = exp(seq(log(lambda_max), log(0.001), length = n_lambda))
# pick the best lambda (chosen via cross-validation)
cv_ridge = cv.glmnet(as.matrix(trainset[ , -response_index]), trainset[ , response_index], alpha = 0, lambda = lambda_seq, nfolds = 5, thresh = 0.0001, standardize = TRUE)
plot(cv_ridge)
# plot(cv_ridge$cvm ~ cv_ridge$lambda, type = "l", ylim = range(cv_ridge$cvm - cv_ridge$cvsd, cv_ridge$cvm + cv_ridge$cvsd))
# abline(v = cv_ridge$lambda.min, col = 2)
text(x = log(cv_ridge$lambda.min), y = 3, paste("lambda:", round(cv_ridge$lambda.min, 3)))
# lines(cv_ridge$cvm + cv_ridge$cvsd ~ cv_ridge$lambda, col = 3)
# lines(cv_ridge$cvm - cv_ridge$cvsd ~ cv_ridge$lambda, col = 3)
# fit final model using ridge, get its sum of squared residuals and multiple R-squared
lambda_cv = cv_ridge$lambda.min
model_ridge = glmnet(as.matrix(trainset[ , -response_index]), trainset[ , response_index], alpha = 0, lambda = lambda_cv, standardize = TRUE)
Yhat_ridge = predict(model_ridge, as.matrix(trainset[ , -response_index]))
rsq_ridge = cor(Ytrain, Yhat_ridge)^2
rss_ridge = sum((Ytrain - Yhat_ridge)^2)
Yhat_test_ridge = predict(model_ridge, as.matrix(testset[ , -response_index]))
rss_test_ridge = sum((Ytest - Yhat_test_ridge)^2)

thresh = 0.2
which_ridge_geqthresh = which(abs(model_ridge$beta) >= thresh) # which |beta|<0.25

mean_beta_full = rep(0, numSat)
beta_true = beta_useful
indices_true = ind_useful
indices0 = which_ridge_geqthresh
indices1 = 1:numSat
mean_beta0 = rep(0, length(indices0))
mean_beta1 = rep(0, length(indices1))

var_e = 1
sigmasq01 = 1

var_beta0 = diag(rep(sigmasq01, length(mean_beta0)))
var_beta1 = diag(rep(sigmasq01, length(mean_beta1)))
numCandidates = 10^5
xmin = 0
xmax = 1
p = length(beta_true)
k = 4 * p
initD = as.matrix(trainset[-response_index])
inity = as.vector(trainset[response_index])
numSeq = 3
N_seq = 10
alpha_seq = 1
buffer_seq = 0
wasserstein0 = 1
genCandidates = 1
seed = 12
