# ridge regression
# https://www.datacamp.com/community/tutorials/tutorial-ridge-lasso-elastic-net

# Load libraries, get data & set seed for reproducibility ---------------------
set.seed(123)    # seef for reproducibility
library(glmnet)  # for ridge regression
library(dplyr)   # for data cleaning
library(psych)   # for function tr() to compute trace of a matrix

data("mtcars")
# Center y, X will be standardized in the modelling function
y <- mtcars %>% select(mpg) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()
X <- mtcars %>% select(-mpg) %>% as.matrix()


# Perform 10-fold cross-validation to select lambda ---------------------------
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
# Setting alpha = 0 implements ridge regression
ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10)
# Plot cross-validation results
plot(ridge_cv)

# Best cross-validated lambda
lambda_cv <- ridge_cv$lambda.min
# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(X, y, alpha = 0, lambda = lambda_cv, standardize = TRUE)
y_hat_cv <- predict(model_cv, X)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2


# Use information criteria to select lambda -----------------------------------
X_scaled <- scale(X)
aic <- c()
bic <- c()
for (lambda in seq(lambdas_to_try)) {
  # Run model
  model <- glmnet(X, y, alpha = 0, lambda = lambdas_to_try[lambda], standardize = TRUE)
  # Extract coefficients and residuals (remove first row for the intercept)
  betas <- as.vector((as.matrix(coef(model))[-1, ]))
  resid <- y - (X_scaled %*% betas)
  # Compute hat-matrix and degrees of freedom
  ld <- lambdas_to_try[lambda] * diag(ncol(X_scaled))
  H <- X_scaled %*% solve(t(X_scaled) %*% X_scaled + ld) %*% t(X_scaled)
  df <- tr(H)
  # Compute information criteria
  aic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df
  bic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df * log(nrow(X_scaled))
}

# Plot information criteria against tried values of lambdas
plot(log(lambdas_to_try), aic, col = "orange", type = "l",
     ylim = c(190, 260), ylab = "Information Criterion")
lines(log(lambdas_to_try), bic, col = "skyblue3")
legend("bottomright", lwd = 1, col = c("orange", "skyblue3"), legend = c("AIC", "BIC"))





####################################################################################
# TOY DATASET ######################################################################
####################################################################################

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
p_useful = 5
p_useless = 100 # not including beta0
p_sat = p_useful + p_useless
beta_useful = runif(n = p_useful, min = 1, max = 2) * sample(c(-1, 1), p_useful, replace = TRUE)
beta_useless = runif(n = p_useless, min = 0, max = 1e-2) * sample(c(-1, 1), p_useless, replace = TRUE)
ind_beta_useful = sort(sample(1:p_sat, p_useful, replace = FALSE))
ind_beta_useless = (1:p_sat)[-ind_beta_useful]
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

# train and test set for cv?

# fit saturated linear model
lmsat = lm(Y ~ X)
rsq_lmsat = invisible(cor(Y, predict(lmsat, data.frame(X)))^2)

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
text(x = cv_lasso$lambda.min, y = 8, paste("lambda:", round(cv_lasso$lambda.min, 3)))
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

## models

library(knitr)
rsqs = round(c(rsq_lmsat, rsq_lmdesired, rsq_lmLASSO), 3)
bics = round(c(BIC(lmsat), BIC(lmdesired), BIC(lmLASSO)), 3)
incl_useful = c(TRUE, TRUE, all(ind_beta_useful %in% lasso_variables))
num_useless = c(p_useless, 0, sum(lasso_variables %in% ind_beta_useless))
p_tot = c(p_sat, p_useful, length(lasso_variables))
df = data.frame(rsqs, bics, incl_useful, num_useless, p_tot, row.names = c("saturated", "sparse", "lasso"))
kable(df)




# https://stats.stackexchange.com/questions/108529/comparing-ols-ridge-and-lasso
