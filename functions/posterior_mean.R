# require("construct_design_matrix.R")
# require("posterior_variance.R")


postmean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       hypothesis_model_type, diagPrior = TRUE){
  X = constructDesignX(D, N, hypothesis_model_type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, hypothesis_model_type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}