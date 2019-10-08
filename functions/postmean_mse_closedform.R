# require("construct_design_matrix.R")
# require("posterior_variance.R")

getClosedEBn = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
  # 1. calculate the variance of the posterior mean
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  # var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% beta_prior_var %*% XtX %*% Sigma_B
  var_postmean_term2 = 0
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1 + var_postmean_term2)
  # 2. calculate expectation of posterior mean
  # expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  expect_postmean = (1/var_e) * Sigma_B %*% XtX %*% matrix(true_beta) + Sigma_B %*% solve(beta_prior_var) %*% matrix(beta_prior_mean)
  return(list("EBn" = expect_postmean, "VarBn" = var_postmean))
}


getClosedMSE = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
  # 1. calculate the variance of the posterior mean
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  # var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% beta_prior_var %*% XtX %*% Sigma_B
  var_postmean_term2 = 0
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1 + var_postmean_term2)
  # 2. calculate expectation of posterior mean
  # expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  expect_postmean = (1/var_e) * Sigma_B %*% XtX %*% matrix(true_beta) + Sigma_B %*% solve(beta_prior_var) %*% matrix(beta_prior_mean)
  # 3. MSE is a function of variance and expectation of posterior mean, as well as the true mean
  biassq_postmean = expect_postmean^2 - 2 * true_beta * expect_postmean + true_beta^2
  MSE_postmean = var_postmean + biassq_postmean
  return(list("var_term" = var_postmean, "biassq_term" = biassq_postmean, "MSE_postmean" = as.vector(MSE_postmean)))
}

# getClosedMSEold = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
#   X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
#   Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
#   # 1. calculate the variance of the posterior mean
#   XtX = crossprod(X)
#   var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
#   var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% beta_prior_var %*% XtX %*% Sigma_B
#   # grab the variance terms in the variance-covariance matrix
#   var_postmean = diag(var_postmean_term1 + var_postmean_term2)
#   # 2. calculate expectation of posterior mean
#   expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
#   # 3. MSE is a function of variance and expectation of posterior mean, as well as the true mean
#   biassq_postmean = expect_postmean^2 - 2 * true_beta * expect_postmean + true_beta^2
#   MSE_postmean = var_postmean + biassq_postmean
#   return(list("var_term" = var_postmean, "biassq_term" = biassq_postmean, "MSE_postmean" = as.vector(MSE_postmean)))
# }