# require("construct_design_matrix.R")

##########
### 2D ###
##########

model_evidence = function(Y, D, N, mean_beta, var_mean, var_e, type){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  if(N != length(Y)) stop("N is not the same length as Y")
  X = constructDesignX(D, N, type)
  marginaly_mean = X %*% mean_beta
  if(dim(X)[1] > 1){
    marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  } else{
    marginaly_var = var_e + (X %*% var_mean %*% t(X))
  }
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}

calcExpPostProbH_2d_old = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
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


calcExpPostProbH_2d = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                   numSims = 100, type = NULL, log_space = TRUE, seed = 123){
  set.seed(seed)
  if(log_space == FALSE){
    # --- Y simulated from H0 --- #
    simY_H0 = simulateY(D, N, mean_beta0, var_mean0, var_e, numSims, type = type[1], seed)
    simEvidenceH0_YH0_elts = matrix(rep(NA, numSims*N), N, numSims)
    simEvidenceH1_YH0_elts = matrix(rep(NA, numSims*N), N, numSims)
    simEvidenceH0_YH0 = rep(NA, numSims)
    simEvidenceH1_YH0 = rep(NA, numSims)
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      for(i in 1:N){
        Y = simY_H0[i, j]
        simEvidenceH0_YH0_elts[i, j] = model_evidence(Y, D[i, , drop = FALSE], 1, mean_beta0, 
                                                    var_mean0, var_e, type = type[1])
        simEvidenceH1_YH0_elts[i, j] = model_evidence(Y, D[i, , drop = FALSE], 1, mean_beta1, 
                                                    var_mean1, var_e, type = type[2])
      }
    }
    simEvidenceH0_YH0 = apply(simEvidenceH0_YH0_elts, 2, prod)
    simEvidenceH1_YH0 = apply(simEvidenceH1_YH0_elts, 2, prod)
    
    # calculate posterior probabilities of models
    simPostH0 = simEvidenceH0_YH0 / (simEvidenceH0_YH0 + simEvidenceH1_YH0)
    simPostH1 = simEvidenceH1_YH0 / (simEvidenceH1_YH0 + simEvidenceH1_YH0)
    # calculate bayes factor
    simBF01 = simPostH0 / simPostH1
    
    expected_postH0_YH0 = mean(simPostH0)
    expected_postH1_YH0 = mean(simPostH1)
    expected_BF01_YH0 = mean(simBF01)
    
    # --- Y simulated from H1 --- #
    simY_H1 = simulateY(D, N, mean_beta1, var_mean1, var_e, numSims, type = type[2], seed)
    simEvidenceH0_YH1_elts = matrix(rep(NA, numSims*N), N, numSims)
    simEvidenceH1_YH1_elts = matrix(rep(NA, numSims*N), N, numSims)
    simEvidenceH0_YH1 = rep(NA, numSims)
    simEvidenceH1_YH1 = rep(NA, numSims)
    simPostH0 = rep(NA, numSims)
    simPostH1 = rep(NA, numSims)
    simBF01 = rep(NA, numSims)
    for(j in 1:numSims){
      for(i in 1:N){
        Y = simY_H1[i, j]
        simEvidenceH0_YH1_elts[i, j] = model_evidence(Y, D[i, , drop = FALSE], 1, mean_beta0, 
                                                      var_mean0, var_e, type = type[1])
        simEvidenceH1_YH1_elts[i, j] = model_evidence(Y, D[i, , drop = FALSE], 1, mean_beta1, 
                                                      var_mean1, var_e, type = type[2])
      }
    }
    simEvidenceH0_YH1 = apply(simEvidenceH0_YH1_elts, 2, prod)
    simEvidenceH1_YH1 = apply(simEvidenceH1_YH1_elts, 2, prod)
    
    # calculate posterior probabilities of models
    simPostH0 = simEvidenceH0_YH1 / (simEvidenceH0_YH1 + simEvidenceH1_YH1)
    simPostH1 = simEvidenceH1_YH1 / (simEvidenceH1_YH1 + simEvidenceH1_YH1)
    # calculate bayes factor
    simBF01 = simPostH0 / simPostH1
    
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

