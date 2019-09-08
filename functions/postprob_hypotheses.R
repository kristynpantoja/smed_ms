model_evidence = function(Y, D, N, mean_beta, var_mean, var_e, type){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  N = length(Y)
  marginaly_mean = X %*% mean_beta
  marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

simulateY = function(D, N, mean_beta, var_mean, var_e, numSims, type = NULL, seed = NULL){
  if(is.null(seed)) set.seed(123)
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    beta = t(rmvnorm(n = 1, mean = mean_beta, sigma = var_mean))
    for(i in 1:N){
      Y[i, j] = rnorm(n = 1, mean = X[i, ] %*% beta, sd = sqrt(var_e))
    }
  }
  return(Y)
}

calcExpPostProbH_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                               numSims = 100, type = NULL, log_space = TRUE, seed = 123){
  set.seed(seed)
  if(log_space == FALSE){
    # --- Y simulated from H0 --- #
    simY = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1], seed)
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
      simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      simBF01[j] = simPostH0[j] / simPostH1[j]
    }
    expected_postH0_YH0 = mean(simPostH0)
    expected_postH1_YH0 = mean(simPostH1)
    expected_BF01_YH0 = mean(simBF01)
    
    # --- Y simulated from H1 --- #
    simY = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
      simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      simBF01[j] = simPostH0[j] / simPostH1[j]
    }
    expected_postH0_YH1 = mean(simPostH0)
    expected_postH1_YH1 = mean(simPostH1)
    expected_BF01_YH1 = mean(simBF01)
  } else{
    # --- Y simulated from H0 --- #
    simY = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
    logsimPostH0 = rep(NA, numSims)
    logsimPostH1 = rep(NA, numSims)
    logsimBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      logsimPostH0[j] = log(simEvidenceH0) - log(simEvidenceH0 + simEvidenceH1)
      logsimPostH1[j] = log(simEvidenceH1) - log(simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      logsimBF01[j] = logsimPostH0[j] - logsimPostH1[j]
    }
    expected_postH0_YH0 = (1 / numSims) * exp(logSumExp(logsimPostH0))
    expected_postH1_YH0 = (1 / numSims) * exp(logSumExp(logsimPostH1))
    expected_BF01_YH0 = (1 / numSims) * exp(logSumExp(logsimBF01))
    
    # --- Y simulated from H1 --- #
    simY = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
    logsimPostH0 = rep(NA, numSims)
    logsimPostH1 = rep(NA, numSims)
    logsimBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      Y = simY[, j]
      # get model evidences
      simEvidenceH0 = model_evidence(Y, D, N, mean_beta0, var_mean0, var_e, type = type[1])
      simEvidenceH1 = model_evidence(Y, D, N, mean_beta1, var_mean1, var_e, type = type[2])
      # calculate posterior probabilities of models
      logsimPostH0[j] = log(simEvidenceH0) - log(simEvidenceH0 + simEvidenceH1)
      logsimPostH1[j] = log(simEvidenceH1) - log(simEvidenceH0 + simEvidenceH1)
      # calculate bayes factor
      logsimBF01[j] = logsimPostH0[j] - logsimPostH1[j]
    }
    expected_postH0_YH1 = (1 / numSims) * exp(logSumExp(logsimPostH0))
    expected_postH1_YH1 = (1 / numSims) * exp(logSumExp(logsimPostH1))
    expected_BF01_YH1 = (1 / numSims) * exp(logSumExp(logsimBF01))
  }
  return(c("expected_postH0_YH0" = expected_postH0_YH0, "expected_postH1_YH0" = expected_postH1_YH0,
           "expected_BF01_YH0" = expected_BF01_YH0, "expected_postH0_YH1" = expected_postH0_YH1,
           "expected_postH1_YH1" = expected_postH1_YH1, "expected_BF01_YH1" = expected_BF01_YH1))
}

