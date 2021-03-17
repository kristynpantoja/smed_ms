# calculate mean-squared error of beta-hat, where beta-hat is the posterior mean
#   of beta
getMSEBeta = function(
  D, N, beta.true, beta.mean, beta.var, error.var, type, indices = NULL, 
  diagPrior = TRUE
  ){
  # design matrix, depends on the model (type)
  X = constructDesignX(D, N, type)
  if(!is.null(type)){
    # posterior variance of beta
    Sigma_B = postvar(D, N, error.var, beta.var, type, diagPrior)
    XtX = crossprod(X)
  } else{
    if(is.null(indices)) stop("if type is NULL, must provide indices")
    # posterior variance
    Sigma_B = postvar(D[ , indices], N, error.var, beta.var, type, diagPrior)
    XtX = crossprod(X[ , indices])
  }
  # 1. calculate the variance of the posterior mean of beta
  var_postmean_term1 = (1/error.var) * Sigma_B %*% XtX %*% Sigma_B
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1)
  # 2. calculate expectation of posterior mean
  expect_postmean = (1/error.var) * Sigma_B %*% XtX %*% matrix(beta.true) + 
    Sigma_B %*% solve(beta.var) %*% matrix(beta.mean)
  # 3. MSE is a function of variance and expectation of posterior mean, 
  #   as well as the true mean
  biassq_postmean = expect_postmean^2 - 2 * beta.true * expect_postmean + beta.true^2
  MSE_postmean = var_postmean + biassq_postmean
  return(list(
    "var_term" = var_postmean, 
    "biassq_term" = biassq_postmean, 
    "MSE_postmean" = as.vector(MSE_postmean)))
}

# getClosedEBn = function(
#   D, N, beta.true, beta.mean, beta.var, error.var, type, 
#   diagPrior = TRUE
#   ){
#   X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
#   # 1. calculate the variance of the posterior mean
#   Sigma_B = postvar(D, N, error.var, beta.var, type, diagPrior)
#   XtX = crossprod(X)
#   var_postmean_term1 = (1/error.var) * Sigma_B %*% XtX %*% Sigma_B
#   # grab the variance terms in the variance-covariance matrix
#   var_postmean = diag(var_postmean_term1)
#   # 2. calculate expectation of posterior mean
#   expect_postmean = (1/error.var) * Sigma_B %*% XtX %*% matrix(beta.true) + 
#     Sigma_B %*% solve(beta.var) %*% matrix(beta.mean)
#   return(list(
#     "EBn" = expect_postmean, 
#     "VarBn" = var_postmean))
# }