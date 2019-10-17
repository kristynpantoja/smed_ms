# require("construct_design_matrix.R")
# require("simulate_y.R")

##########
### 2D ###
##########

model_evidence = function(Y, D, N, beta_prior_mean, beta_prior_var, var_e, hypothesis_model_type){
  # Y is a vector of outputs (or scalar, for one output)
  # X is a matrix of inputs (or vector, for one input)
  # beta_prior_var is a matrix
  # var_e is a scalar
  if(N != length(Y)) stop("N is not the same length as Y")
  X = constructDesignX(D, N, hypothesis_model_type)
  # get mean and variance of marginal density of y, which is N-dim multivariate normal pdf
  marginaly_mean = X %*% beta_prior_mean
  if(dim(X)[1] > 1){ # if X is a matrix of inputs
    marginaly_var = diag(rep(var_e, N)) + (X %*% beta_prior_var %*% t(X))
  } else{ # if X is a vector for one input
    marginaly_var = var_e + (X %*% beta_prior_var %*% t(X))
  }
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

calcExpPostProbH = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
                            beta_prior_var0, beta_prior_var1, var_e,
                            numSims = 100, true_model_type = NULL, H01_model_types = NULL,
                            seed = 123){
  # for now we assume that the type of model for true_beta contains the types for both hypotheses
  simY = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
  simPostH0 = rep(NA, numSims)
  simPostH1 = rep(NA, numSims)
  simBF01 = rep(NA, numSims)
  for(j in 1:numSims){
    Y = simY[ , j]
    # get model evidences for each hypothesized model
    simEvidenceH0 = model_evidence(Y, D, N, beta_prior_mean0, beta_prior_var0, var_e, H01_model_types[1])
    simEvidenceH1 = model_evidence(Y, D, N, beta_prior_mean1, beta_prior_var1, var_e, H01_model_types[2])
    # calculate posterior probabilities of each hypothesis
    simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
    simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
    # calculate bayes factor
    simBF01[j] = simPostH0[j] / simPostH1[j]
  }
  expected_postH0 = mean(simPostH0)
  expected_postH1 = mean(simPostH1)
  expected_BF01 = mean(simBF01)
  
  return(c("expected_postH0" = expected_postH0, "expected_postH1" = expected_postH1,
           "expected_BF01" = expected_BF01))
}

calcExpPostProbH_data = function(y, D, N, beta_prior_mean0, beta_prior_var0, 
                                 beta_prior_mean1, beta_prior_var1, var_e, model_types){
  # get model evidence
  evidence0 = model_evidence(y, D, N, beta_prior_mean0, beta_prior_var0, var_e, model_types[1])
  evidence1 = model_evidence(y, D, N, beta_prior_mean1, beta_prior_var1, var_e, model_types[2])
  # get each hypotheses' posterior probability
  # calculate posterior probabilities of each hypothesis
  postprob0 = evidence0 / (evidence0 + evidence1)
  postprob1 = evidence1 / (evidence0 + evidence1)
  BF01 = evidence0 / evidence1
  return(c("postprob0" = postprob0, "postprob1" = postprob1, "BF01" = BF01))
}





# mean_beta0 = beta_prior_mean0
# mean_beta1 = beta_prior_mean1
# var_mean0 = beta_prior_var0
# var_mean1 = beta_prior_var1
# calcExpPostProbH_2d_old = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
#                                numSims = 100, type = NULL, seed = 123){
#   set.seed(seed)
#   # --- Y simulated from H0 --- #
#   simY = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1], seed)
#   simPostH0 = rep(NA, numSims)
#   simPostH1 = rep(NA, numSims)
#   simBF01 = rep(NA, numSims)
#   for(j in 1:numSims){
#     Y = simY[, j]
#     # get model evidences for each hypothesized model
#     simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
#     simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
#     # calculate posterior probabilities of each hypothesis
#     simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
#     simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
#     # calculate bayes factor
#     simBF01[j] = simPostH0[j] / simPostH1[j]
#   }
#   expected_postH0_YH0 = mean(simPostH0)
#   expected_postH1_YH0 = mean(simPostH1)
#   expected_BF01_YH0 = mean(simBF01)
#   
#   # --- Y simulated from H1 --- #
#   simY = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
#   simPostH0 = rep(NA, numSims)
#   simPostH1 = rep(NA, numSims)
#   simBF01 = rep(NA, numSims)
#   for(j in 1:numSims){
#     Y = simY[, j]
#     # get model evidences
#     simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
#     simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
#     # calculate posterior probabilities of models
#     simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
#     simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
#     # calculate bayes factor
#     simBF01[j] = simPostH0[j] / simPostH1[j]
#   }
#   expected_postH0_YH1 = mean(simPostH0)
#   expected_postH1_YH1 = mean(simPostH1)
#   expected_BF01_YH1 = mean(simBF01)
#   
#   return(c("expected_postH0_YH0" = expected_postH0_YH0, "expected_postH1_YH0" = expected_postH1_YH0,
#            "expected_BF01_YH0" = expected_BF01_YH0, "expected_postH0_YH1" = expected_postH0_YH1,
#            "expected_postH1_YH1" = expected_postH1_YH1, "expected_BF01_YH1" = expected_BF01_YH1))
# }
# 
# 

