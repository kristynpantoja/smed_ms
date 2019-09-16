# require("construct_design_matrix.R")
# require("simulate_y.R")

# first compute estimator, posterior mean
getPostMean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       hypothesis_model_type, diagPrior = TRUE){
  X = constructDesignX(D, N, hypothesis_model_type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, hypothesis_model_type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}

getEmpMSE = function(postmean_beta, true_beta){
  return(sum((postmean_beta - true_beta)^2))
}

calcExpPostMeanMSE = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
                              beta_prior_var0, beta_prior_var1, var_e,
                              numSims = 100, true_model_type = NULL, H01_model_types = NULL, 
                              seed = 123, diagPrior = TRUE){
  set.seed(seed)
  Ysims = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
  # number of parameters for beta in hypothesized models 
  p0 = length(beta_prior_mean0)
  p1 = length(beta_prior_mean1)
  p_true = length(true_beta)
  # calculating posterior means from each hypothesis' prior on beta
  if(true_model_type == 5 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
    # posterior mean, given H0 prior on beta
    p0_indices_wrt_H1 = c(1, 2, 4) # indices of parameters in null model, i.e. intercept and linear terms
    post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
                                                                 var_e, type[1], diagPrior))
    empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getEmpMSE(x, true_beta[p0_indices_wrt_H1])) # Only MSE if assume H0 for estimation
    expEmpiricalMSEH0 = mean(empMSEsH0)
    # posterior mean, given H1 prior on beta
    post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1, beta_prior_var1, 
                                                                 var_e, type[2], diagPrior))
    empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getEmpMSE(x, true_beta))
    expEmpiricalMSEH1 = mean(empMSEsH1)
  } else if(true_model_type == 4 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
    # posterior mean, given H0 prior on beta
    p1_indices_wrt_H0 = c(1, 2, 4) # indices of parameters in alt model, i.e. intercept and linear terms
    post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
                                                                 var_e, type[1], diagPrior))
    empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getEmpMSE(x, true_beta)) # Only MSE if assume H0 for estimation
    expEmpiricalMSEH0 = apply(empMSEsH0, 1, mean)
    # posterior mean, given H1 prior on beta
    post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1[p1_indices_wrt_H0], 
                                                                 diag(diag(beta_prior_var1)[p1_indices_wrt_H0]), 
                                                                 var_e, type[1], diagPrior))
    empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getEmpMSE(x, true_beta))
    expEmpiricalMSEH1 = apply(empMSEsH1, 1, mean)
  } else{
    expEmpiricalMSEH0 = rep(NA, p0)
    expEmpiricalMSEH1 = rep(NA, p1)
  }
  # # put the results together
  # expEmpiricalMSEs = matrix(NA, p_true, 2)
  # if(true_model_type == 5 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
  #   expEmpiricalMSEs[p0_indices_wrt_H1, 1] = expEmpiricalMSEH0
  #   expEmpiricalMSEs[1:p1, 2] = expEmpiricalMSEH1
  #   colnames(expEmpiricalMSEs) = c("expEmpiricalMSEH0", "expEmpiricalMSEH1")
  #   rownames(expEmpiricalMSEs) = paste("E[MSE(B", 0:(dim(expEmpiricalMSEs)[1] - 1), ")]", sep = "")
  # } else if(true_model_type == 4 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
  #   expEmpiricalMSEs[1:p0, 1] = expEmpiricalMSEH0
  #   expEmpiricalMSEs[1:p0, 2] = expEmpiricalMSEH1
  #   colnames(expEmpiricalMSEs) = c("expEmpiricalMSEH0", "expEmpiricalMSEH1")
  #   rownames(expEmpiricalMSEs) = paste("E[MSE(B", 0:(dim(expEmpiricalMSEs)[1] - 1), ")]", sep = "")
  # }
  return(c("expEmpiricalMSEH0" = expEmpiricalMSEH0, "expEmpiricalMSEH1" = expEmpiricalMSEH1))
}


# calcExpPostMeanMSE_old = function(D, N, true_beta, numSims, mean_beta, var_e, var_mean, type, diagPrior = TRUE){
#   Ysims = simulateY(D, N, mean_beta, var_mean, var_e, numSims, type = type)
#   post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, true_beta, var_e, var_mean, type, diagPrior = diagPrior))
#   empMSEs = apply(post_means, 2, FUN = function(x) getEmpMSE(x, true_beta))
#   expEpiricalMSE = apply(empMSEs, 1, mean)
#   return(expEpiricalMSE)
# }
