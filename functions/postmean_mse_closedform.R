# require("construct_design_matrix.R")

##########
### 2D ###
##########

# closed-form of MSE of posterior mean = variance of posterior mean if true_beta = mean_beta (since no bias)
getClosedMSE = function(D, N, true_beta, mean_beta, var_e, var_mean, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  Sigma_B = postvar(D, N, var_e, var_mean, type, diagPrior) # posterior variance
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% var_mean %*% XtX %*% Sigma_B
  var_postmean = var_postmean_term1 + var_postmean_term2
  expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(var_mean)) %*% matrix(mean_beta)
  MSE_postmean = diag(var_postmean) + expect_postmean^2 - 2 * true_beta * expect_postmean + true_beta^2
  return(as.vector(MSE_postmean))
}
