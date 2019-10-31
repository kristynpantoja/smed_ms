# require("construct_design_matrix.R")
# require("posterior_variance.R")

getClosedMSEyhat = function(pt, D, N, true_beta, true_type, beta_prior_mean, beta_prior_var, var_e, prior_type, diagPrior = TRUE){
  # note true_type must match prior_type
  if(true_type != prior_type) warning("true_type is not the same as prior_type: will lead to non-conforming matrices in multiplication")
  # variance term
  x.temp = matrix(constructDesignX(pt, 1, prior_type))
  X.temp = constructDesignX(D, N, prior_type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, prior_type) # posterior variance
  XtX = crossprod(X.temp)
  varBn = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  variance_term = t(x.temp) %*% varBn %*% x.temp
  
  # bias-squared term
  # calculate expectation of posterior mean
  # expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  expect_postmean = (1/var_e) * Sigma_B %*% XtX %*% matrix(true_beta) + Sigma_B %*% solve(beta_prior_var) %*% matrix(beta_prior_mean)
  x.true = matrix(constructDesignX(pt, 1, true_type))
  biassq_term = (t(x.temp) %*% expect_postmean - t(x.true) %*% true_beta)^2
  
  # MSE of y-hat
  MSEyhat = variance_term + biassq_term
  return(list("var_term" = variance_term, "biassq_term" = biassq_term, "MSEyhat" = as.vector(MSEyhat)))
}
getClosedMSEyhat_seq = function(x_seq, D, N, true_beta, true_type, beta_prior_mean, beta_prior_var, var_e, prior_type, diagPrior = TRUE){
  mseyhat_seq = rep(NA, length(x_seq))
  mseyhat_var_seq = rep(NA, length(x_seq))
  mseyhat_biassq_seq = rep(NA, length(x_seq))
  for(i in 1:length(x_seq)){
    mseyhat = getClosedMSEyhat(x_seq[i], D, N, true_beta, true_type, beta_prior_mean, beta_prior_var, var_e, prior_type, diagPrior)
    mseyhat_seq[i] = mseyhat$MSEyhat
    mseyhat_var_seq[i] = mseyhat$var_term
    mseyhat_biassq_seq[i] = mseyhat$biassq_term
  }
  return(list("MSEyhat" = mseyhat_seq, "var_term" = mseyhat_var_seq, "biassq_term" = mseyhat_biassq_seq))
}
