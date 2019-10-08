# require("construct_design_matrix.R")
# require("posterior_variance.R")


postmean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       hypothesis_model_type, diagPrior = TRUE){
  X = constructDesignX(D, N, hypothesis_model_type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, hypothesis_model_type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}

calcExpPostMean = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
                           beta_prior_var0, beta_prior_var1, var_e,
                           numSims = 100, true_model_type = NULL, H01_model_types = NULL,
                           diagPrior = TRUE, seed = 123){
  simY = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
  Bns0 = apply(simY, 2, FUN = function(y) postmean(y, D, N, beta_prior_mean0, beta_prior_var0,
                                                        var_e, H01_model_types[1], diagPrior))
  expBns0 = apply(Bns0, 1, mean)
  var0 = apply(Bns0, 1, var)
  Bns1 = apply(simY, 2, FUN = function(y) postmean(y, D, N, beta_prior_mean1, beta_prior_var1,
                                                   var_e, H01_model_types[2], diagPrior))
  expBns1 = apply(Bns1, 1, mean)
  var1 = apply(Bns1, 1, var)
  return(list("expBns0" = expBns0, "expBns1" = expBns1, "var0" = var0, "var1" = var1))
}
