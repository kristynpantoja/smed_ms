# require("construct_design_matrix.R")
# require("simulate_y.R")

model_evidence = function(Y, D, N, beta_prior_mean, beta_prior_var, var_e, hypothesis_model_type, indices = NULL){
  # Y is a vector of outputs (or scalar, for one output)
  # X is a matrix of inputs (or vector, for one input)
  # beta_prior_var is a matrix
  # var_e is a scalar
  if(N != length(Y)) stop("N is not the same length as Y")
  if(!is.null(hypothesis_model_type)){
    X = constructDesignX(D, N, hypothesis_model_type)
  } else{
    X = constructDesignX(D, N, hypothesis_model_type)[ , indices]
  }
  # get mean and variance of marginal density of y, which is N-dim multivariate normal pdf
  marginaly_mean = X %*% beta_prior_mean
  if(dim(X)[1] > 1){ # if X is a matrix of inputs
    marginaly_var = diag(rep(var_e, N)) + (X %*% beta_prior_var %*% t(X))
  } else{ # if X is a vector for one input
    warning("X is a vector, not a matrix - is that what you expected?")
    marginaly_var = var_e + (X %*% beta_prior_var %*% t(X))
  }
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

calcExpPostProbH = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
                            beta_prior_var0, beta_prior_var1, var_e,
                            numSims = 100, true_model_type = NULL, H01_model_types = NULL,
                            seed = NULL){
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
    simBF01[j] = simEvidenceH0 / simEvidenceH1
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



## for M hypotheses

# for non-sequential designs, where data needs to be generated in simulations to estimate
calcEPPH = function(D, N, true_beta, true_model_type, models, var_e, numSims = 100, true_indices = NULL, seed = NULL){
  if(!is.null(true_model_type)){
    simY = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
  } else{
    simY = simulateYvs(D[ , true_indices], N, true_beta, var_e, numSims, seed)
  }
  # "models" is a list of lists, where each element of list "models" describes a model
  #   each element/model is a list containing "beta_prior_mean", "beta_prior_var", "model_type"
  #   in that order
  #   if the model type is for variable selection, i.e. model_type == NULL, there will be a fourth element
  #   in the model's list: indices, which corresponds to its hypothesized beta_prior_mean
  # for now we assume that the type of model for true_beta contains the types for both hypotheses
  # model_evidences = matrix(NA, length(models), numSims)
  model_postprobs = matrix(NA, length(models), numSims)
  for(j in 1:numSims){
    Y = simY[ , j]
    model_postprobs[ , j] = calcEPPHdata(Y, D, N, models, var_e)
  }
  exp_postprobs = apply(model_postprobs, 1, mean)
  
  return(c("exp_postprobs" = exp_postprobs))
}

# for designs with data
calcEPPHdata = function(y, D, N, models, var_e){
  model_evidences = rep(NA, length(models))
  model_postprobs = rep(NA, length(models))
  # get model evidence
  for(m in 1:length(models)){
    model = models[[m]]
    if(length(model) == 3) model[[4]] = NULL
    model_evidences[m] = model_evidence(y, D, N, model[[1]], model[[2]], var_e, model[[3]], model[[4]])
  }
  # get each hypotheses' posterior probability
  denom = sum(model_evidences)
  model_postprobs = model_evidences / denom
  return(c("model_postprobs" = model_postprobs))
}

calcEPPHseqdata = function(y, D, models, var_e, initN, numSeq, N_seq){
  postprobs = matrix(NA, length(models), numSeq)
  for(i in 1:numSeq){
    changing_postprobs = calcEPPHdata(y[1:(initN + N_seq * i)],
                                      D[1:(initN + N_seq * i), ], 
                                      N = initN + N_seq * i, models, sigmasq)
    postprobs[ , i] = changing_postprobs
  }
  return(postprobs)
}
