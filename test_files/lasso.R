# LASSO

# Load the riboflavin data
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression
#dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
#?riboflavin # this gives you more information on the dataset

X = riboflavin$x
Y = riboflavin$y

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

# get lambda_seq
n = dim(X)[1]
p = dim(X)[2]
XYstd = standardizeXY(X, Y)
Xtilde = XYstd$Xtilde
Ytilde = XYstd$Ytilde
n_lambda = 30
lambda_max = max(abs(crossprod(Xtilde, Ytilde) / n))
lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

library(glmnet)

cv_lasso = cv.glmnet(Xtilde, Ytilde, alpha = 1, lambda = lambda_seq, nfolds = 5, thresh = 0.0001, standardize = TRUE)
plot(cv_lasso$cvm ~ lambda_seq, type = "l")
abline(v = cv_lasso$lambda.min, col = 2)
abline(v = cv_lasso$lambda.1se, col = 2)
lines(cv_lasso$cvm + cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)
lines(cv_lasso$cvm - cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)

# print(paste("glmnet lambda_min: ", cv_lasso$lambda.min))
# print(paste("glmnet lambda_1se: ", cv_lasso$lambda.1se))

# Best cross-validated lambda
lasso_cv = cv.glmnet(Xtilde, Ytilde, alpha = 1, lambda = lambda_seq, standardize = TRUE, nfolds = 10)
# Fit final model, get its sum of squared residuals and multiple R-squared
lambda_cv = cv_lasso$lambda.min
model_cv = glmnet(Xtilde, Ytilde, alpha = 1, lambda = lambda_cv, standardize = TRUE)
y_hat_cv = predict(model_cv, Xtilde)
ssr_cv = t(Ytilde - y_hat_cv) %*% (Ytilde - y_hat_cv)
rsq_lasso_cv = cor(Ytilde, y_hat_cv)^2

# lmSaturated = lm(Ytilde ~ Xtilde) # saturated model is not valid!
# rsq_lmSaturated = cor(Ytilde, predict(lmSaturated, Xtilde))^2
# lmLASSOvars

##############################################################################################
##############################################################################################
##############################################################################################
# TOY DATA SET #
data(QuickStartExample)
n = dim(x)[1]
p = dim(x)[2]
# standardize and get lambda_seq
xystd = standardizeXY(x, y)
xtilde = xystd$Xtilde
ytilde = xystd$Ytilde
n_lambda = 30
lambda_max = max(abs(crossprod(xtilde, ytilde) / n))
lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

cv_lasso = cv.glmnet(xtilde, ytilde, alpha = 1, lambda = lambda_seq, nfolds = 5, thresh = 0.0001, standardize = TRUE)
plot(cv_lasso$cvm ~ lambda_seq, type = "l")
abline(v = cv_lasso$lambda.min, col = 2)
abline(v = cv_lasso$lambda.1se, col = 2)
lines(cv_lasso$cvm + cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)
lines(cv_lasso$cvm - cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)

# print(paste("glmnet lambda_min: ", cv_lasso$lambda.min))
# print(paste("glmnet lambda_1se: ", cv_lasso$lambda.1se))

# Best cross-validated lambda
lasso_cv = cv.glmnet(xtilde, ytilde, alpha = 1, lambda = lambda_seq, standardize = TRUE, nfolds = 10)
# Fit final model, get its sum of squared residuals and multiple R-squared
lambda_cv = cv_lasso$lambda.min
model_cv = glmnet(xtilde, ytilde, alpha = 1, lambda = lambda_cv, standardize = TRUE)
y_hat_cv = predict(model_cv, xtilde)
ssr_cv = t(ytilde - y_hat_cv) %*% (ytilde - y_hat_cv)
rsq_lasso_cv = cor(ytilde, y_hat_cv)^2

lmSaturated = lm(ytilde ~ xtilde)
rsq_lmSaturated = cor(ytilde, predict(lmSaturated, data.frame(xtilde)))^2

lasso_variables = which(as.matrix(model_cv$beta) != 0)
lasso_xtilde = xtilde[ , lasso_variables]
lmLASSOvars = lm(ytilde ~ lasso_xtilde)
rsq_lmLASSOvars = cor(ytilde, predict(lmLASSOvars, data.frame(lasso_xtilde)))^2

# want to choose x to better be able to distinguish these two models
##############################################################################################
##############################################################################################
##############################################################################################
# TOY DATA SET #

# libraries and helper fns
library(mvtnorm)
library(glmnet)
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
p_useful = 5
p_useless = 100 # not including beta0
p_sat = p_useful + p_useless
beta_useful = runif(n = p_useful, min = 1, max = 2) * sample(c(-1, 1), p_useful, replace = TRUE)
beta_useless = runif(n = p_useless, min = 0, max = 1e-2) * sample(c(-1, 1), p_useless, replace = TRUE)
ind_beta_useful = sort(sample(1:p_sat, p_useful, replace = FALSE))
beta_sat = rep(NA, p_sat)
beta_sat[ind_beta_useful] = sample(beta_useful)
beta_sat[-ind_beta_useful] = sample(beta_useless)
# beta_sat = sample(c(beta_useful, beta_useless), p_sat, replace = FALSE)
plot(x = 1:p_sat, y = beta_sat); abline(h = 0, lty = 3)

# define true model
beta0 = 0
truef = function(x) beta0 + t(beta_sat) %*% x

# generate data
set.seed(123*sample(1:50000, 1))
xmin = -1
xmax = 1
Sigma = matrix(0.5, p_sat, p_sat) + diag(rep(1-0.5, p_sat)) # covariance for design X
n = 100 # sample size
X = rmvnorm(n, mean = rep(0, p_sat), sigma = Sigma) # matrix of covariates
sigma = 1 # noise standard deviation
Y = beta0 + X %*% beta_sat + sigma * rnorm(n) # response from linear model

# train and test set


# fit saturated linear model
lmsat = lm(Y ~ X)
rsq_lmsat = cor(Y, predict(lmsat, data.frame(X)))^2

# fit desired model (no useless variables)
Xdesired = X[ , ind_beta_useful]
lmdesired = lm(Y ~ Xdesired)
rsq_lmdesired = cor(Y, predict(lmdesired, data.frame(Xdesired)))^2


# fit LASSO linear model
# get lambda_seq
n = dim(X)[1]
p = dim(X)[2]
XYstd = standardizeXY(X, Y)
Xtilde = XYstd$Xtilde
Ytilde = XYstd$Ytilde
n_lambda = 30
lambda_max = max(abs(crossprod(Xtilde, Ytilde) / n))
lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
cv_lasso = cv.glmnet(Xtilde, Ytilde, alpha = 1, lambda = lambda_seq, nfolds = 5, thresh = 0.0001, standardize = TRUE)
plot(cv_lasso$cvm ~ lambda_seq, type = "l")
abline(v = cv_lasso$lambda.min, col = 2)
abline(v = cv_lasso$lambda.1se, col = 2)
lines(cv_lasso$cvm + cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)
lines(cv_lasso$cvm - cv_lasso$cvsd ~ cv_lasso$lambda, col = 3)
# get best cross-validated lambda
lasso_cv = cv.glmnet(X, Y, alpha = 1, lambda = lambda_seq, standardize = TRUE, nfolds = 10)
# fit final model using LASSO, get its sum of squared residuals and multiple R-squared
lambda_cv = cv_lasso$lambda.min
model_lasso = glmnet(X, Y, alpha = 1, lambda = lambda_cv, standardize = TRUE)
y_hat_lasso = predict(model_lasso, X)
ssr_lasso = t(Y - y_hat_lasso) %*% (Y - y_hat_lasso)
rsq_lasso = cor(Y, y_hat_lasso)^2
# # backscale (didn't implement bc not necessary for glmnet when standardize = TRUE)
# # Perform back scaling and centering to get original intercept and coefficient vector for each lambda
# lasso_variables = which(as.vector(model_lasso$beta) != 0)
# beta_lasso = diag(XYstd$weights) %*% model_lasso$beta
# beta0_lasso = XYstd$Ymean - crossprod(XYstd$Xmeans, beta)
# fit linear model using LASSO variables
lasso_variables = which(as.vector(model_lasso$beta) != 0)
X_lasso = X[ , lasso_variables]
lmLASSO = lm(Y ~ X_lasso)
rsq_lmLASSO = cor(Y, predict(lmLASSO, data.frame(X_lasso)))^2

# check
# saturated model
rsq_lmsat
BIC(lmsat)
# desired model
rsq_lmdesired
BIC(lmdesired)
# lasso
rsq_lmLASSO
BIC(lmLASSO)
p_sat - length(lasso_variables) # ideally, this would be close to p_useless = 200
all(ind_beta_useful %in% lasso_variables) # that all the useful variables are included
# it's okay if some of the useless ones are in there, I guess. maybe they're not actually useless.

# what points are optimal for distinguishing desired model and saturated model?
# what points are optimal for distinguishing LASSO model and saturated model?
# what regions does LASSO have an easier time choosing the best parameters beta?

# apply mmed and smmed just to see what happens, I guess



# 
# # first check smmed fn
# # libraries
# library(expm)
# library(matrixStats)
# library(scatterplot3d)
# library(knitr)
# library(mvtnorm)
# 
# # source files for evaluations
# 
# # --- sources to generate MEDs --- #
# home = "/home/kristyn/Documents/smed_ms"
# functions_home = paste(home, "/functions", sep="")
# source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
# source(paste(functions_home, "/charge_function_q.R", sep = ""))
# source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
# source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
# source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
# 
# source(paste(functions_home, "/posterior_mean.R", sep = ""))
# source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
# source(paste(functions_home, "/posterior_variance.R", sep = ""))
# source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
# source(paste(functions_home, "/simulate_seqMED.R", sep = ""))
# 
# # --- sources to designs : MSE(Bn), E[P(H1|Y,D)] --- #
# source(paste(functions_home, "/simulate_y.R", sep = ""))
# source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
# source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
# source(paste(functions_home, "/plot_EPH1.R", sep = ""))
# source(paste(functions_home, "/plot_MSE.R", sep = ""))
# source(paste(functions_home, "/plot_posterior_variance.R", sep = ""))
# source(paste(functions_home, "/postpredyhat_mse_closedform.R", sep = ""))
# 
# # CASE 1 #
# 
# mu0 = c(0, 0)
# mu1 = c(0, 0, 0)
# typeT = 3
# betaT = c(-0.2, -0.4, 0.4)
# sigmasq01 = 0.25
# sigmasq = 0.1
# 
# # MED design #
# typeT = 3
# f0 = function(x) mu0[1] + mu0[2] * x
# f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
# xmin = -1
# xmax = 1
# 
# fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
# curve(f0, col = 2, lwd = 5, xlim = c(xmin, xmax), xlab = "x", ylab = "f")
# curve(f1, add = T, col = 5, lty = 2, lwd = 5)
# curve(fT, add = T, col = 1, lty = 3, lwd = 5)
# legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))
# 
# # MED design #
# V0 = diag(rep(sigmasq01,length(mu0)))
# V1 = diag(rep(sigmasq01,length(mu1)))
# # function settings (including and based on prior settings above)
# f0 = function(x) mu0[1] + mu0[2] * x
# f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
# type01 = c(2, 3)
# numCandidates = 10^3
# k = 4
# xmin = -1
# xmax = 1
# p = 1
# N = 15
# initx = rnorm(5)
# inity = f1(initx) + rnorm(5, 0, sqrt(sigmasq))
# 
# numSeq = 3
# N_seq = 5
# alpha_seq = 1
# smmed4 = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
#                          numCandidates, k, xmin, xmax, p, initD = initx, inity = inity,
#                          numSeq = numSeq, N_seq = N_seq, alpha_seq = alpha_seq, seed = 12)
# 
# smmed5 = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
#                          numCandidates, k, xmin, xmax, p, initD = initx, inity = inity,
#                          numSeq = numSeq - 1, N_seq = N_seq, alpha_seq = alpha_seq, seed = 12)
# 
# length(smmed5$y)
# length(smmed5$D)
