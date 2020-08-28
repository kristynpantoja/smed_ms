# require("construct_design_matrix.R")
# require("posterior_variance.R")

getClosedEBn = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
  # 1. calculate the variance of the posterior mean
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1)
  # 2. calculate expectation of posterior mean
  # expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  expect_postmean = (1/var_e) * Sigma_B %*% XtX %*% matrix(true_beta) + Sigma_B %*% solve(beta_prior_var) %*% matrix(beta_prior_mean)
  return(list("EBn" = expect_postmean, "VarBn" = var_postmean))
}


getClosedMSE = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, indices = NULL, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  if(!is.null(type)){
    Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
    # 1. calculate the variance of the posterior mean
    XtX = crossprod(X)
  } else{
    if(is.null(indices)) stop("if type is NULL, must provide indices")
    Sigma_B = postvar(D[ , indices], N, var_e, beta_prior_var, type, diagPrior) # posterior variance
    # 1. calculate the variance of the posterior mean
    XtX = crossprod(X[ , indices])
  }
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1)
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



#####################################################
# Monte Carlo Approximation, approaches closed form #
#####################################################

# require("construct_design_matrix.R")
# require("simulate_y.R")
# require("posterior_variance.R")

# first compute estimator, posterior mean
getPostMean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       hypothesis_model_type, diagPrior = TRUE){
  X = constructDesignX(D, N, hypothesis_model_type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, hypothesis_model_type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}

getDevianceSq = function(postmean_beta, true_beta){
  return((postmean_beta - true_beta)^2)
}

calcExpPostMeanMSE = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, numSims, diagPrior = TRUE, seed = 123){
  # assumes true_beta and beta_prior_mean have the same model type (linear or quadratic)
  set.seed(seed)
  Ysims = simulateY(D, N, true_beta, var_e, numSims, type, seed)
  # calculating posterior mean
  post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean, beta_prior_var, 
                                                             var_e, type, diagPrior))
  # calculate squared deviances for each parameter Bni in Bn = (Bn1, ..., Bnp) from the true Bi in B = (B1, ..., Bp)
  empMSEs = apply(post_means, 2, FUN = function(x) getDevianceSq(x, true_beta))
  # get the mean squared deviances
  expEmpiricalMSE = apply(empMSEs, 1, mean)
  
  # to look at the variance of posterior mean:
  # var_postmeans = apply(post_means, 1, var)
  # and their means, to see if they're centered, i.e. to see if I calculated Bn correctly
  # avg_postmeans = apply(post_means, 1, mean)
  # return(list("expEmpiricalMSE" = expEmpiricalMSE, "var_postmeans" = var_postmeans, "avg_postmeans" = avg_postmeans))
  return("expEmpiricalMSE" = expEmpiricalMSE)
}

# calcExpPostMeanMSE_old = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
#                               beta_prior_var0, beta_prior_var1, var_e,
#                               numSims = 100, true_model_type = NULL, H01_model_types = NULL, 
#                               seed = 123, diagPrior = TRUE){
#   set.seed(seed)
#   Ysims = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
#   # calculating posterior means from each hypothesis' prior on beta
#   if(true_model_type == 5 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
#     # posterior mean, given H0 prior on beta
#     # we assume linear model to calculate MSE, since this model doesn't give estimates of quadratic terms
#     post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
#                                                                  var_e, type[1], diagPrior))
#     empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getDevianceSq(x, true_beta[c(1, 2, 4)])) # Only MSE if assume linear model for estimation
#     # expEmpiricalMSEH0 = mean(empMSEsH0)
#     expEmpiricalMSEH0 = apply(empMSEsH0, 1, mean)
#     # posterior mean, given H1 prior on beta
#     post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1, beta_prior_var1, 
#                                                                  var_e, type[2], diagPrior))
#     empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getDevianceSq(x, true_beta))
#     # expEmpiricalMSEH1 = mean(empMSEsH1)
#     expEmpiricalMSEH1 = apply(empMSEsH1, 1, mean)
#   } else if(true_model_type == 4 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
#     # posterior mean, given H0 prior on beta
#     post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
#                                                                  var_e, type[1], diagPrior))
#     # get squared deviances for each parameter Bi
#     empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getDevianceSq(x, true_beta)) # Only MSE if assume linear model for estimation
#     expEmpiricalMSEH0 = mean(empMSEsH0) # average squared deviances of each Bi
#     # posterior mean, given H1 prior on beta
#     # post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1[c(1, 2, 4)], 
#     #                                                              diag(diag(beta_prior_var1)[c(1, 2, 4)]), 
#     #                                                              var_e, type[1], diagPrior))
#     # we ignore estimates for quadratic terms
#     post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1, beta_prior_var1, 
#                                                                  var_e, type[2], diagPrior))
#     post_meansH1 = post_meansH1[c(1, 2, 4), ]
#     empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getDevianceSq(x, true_beta))
#     expEmpiricalMSEH1 = mean(empMSEsH1)
#   } else{
#     expEmpiricalMSEH0 = rep(NA, length(beta_prior_mean0))
#     expEmpiricalMSEH1 = rep(NA, length(beta_prior_mean1))
#   }
#   return(list("expEmpiricalMSEH0" = expEmpiricalMSEH0, "expEmpiricalMSEH1" = expEmpiricalMSEH1))
# }