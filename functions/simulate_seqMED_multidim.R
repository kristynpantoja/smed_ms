# for multidimensional case... like LASSO or ridge
simulate_seqMED_multidim = function(mean_beta_full, beta_true, indices_true, indices0, indices1, 
                                    mean_beta0 = NULL, mean_beta1 = NULL, var_beta0, var_beta1, 
                                    xmin = -1, xmax = 1, numCandidates = 10^5, k = 4, p = 1,
                                    initD = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
                                    alpha_seq = NULL, buffer_seq = 0, candidates = NULL, 
                                    wasserstein0 = 1, genCandidates = 1, seed = NULL){
  # some checks!!!!!!!!!!!!!!!!!!!!!!
  if(is.null(initD)){
    stop("no initial data! haven't written this case for multiple dimensions yet.")
    # if(is.null(alpha_seq)){ # generate space-filling design for first step
    #   D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
    #                          numCandidates, k, xmin, xmax, p, alpha = 0, buffer_seq[1], 
    #                          genCandidates = 1, initialpt = 1)
    # } else{
    #   D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
    #                          numCandidates, k, xmin, xmax, p, alpha_seq[1], buffer_seq[1], 
    #                          genCandidates = 1, initialpt = 1)
    # }
    # initN = length(D1)
  } else{
    D1 = initD
    initN = dim(D1)[1]
  }
  
  if(numSeq == 1){
    return(D1)
  }
  # otherwise, sequence:
  # some checks first:
  if(length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
  if(length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(!is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
  if(is.null(initD)){
    Nttl = sum(N_seq)
  } else{
    Nttl = sum(initN, N_seq)
  }
  
  if(!is.null(seed)) set.seed(seed)
  # get y1
  if(is.null(inity)){
    y1 = as.vector(simulateY(D1, N_seq[1], true_beta, sigmasq, 1, true_type))
  } else{
    y1 = inity
  }
  y = y1
  D = D1
  
  # why do I do this? why only take the diagonal of postvar?
  current_postvar0 = diag(postvar(D[ , indices0], dim(D[ , indices0])[1], var_e, var_beta0))
  current_postmean0 = postmean(y, D[ , indices0], dim(D[ , indices0])[1], mean_beta0, var_beta0, var_e, type[1])
  current_postvar1 = diag(postvar(D[ , indices1], dim(D[ , indices1])[1], var_e, var_beta1, 1))
  current_postmean1 = postmean(y, D[ , indices1], dim(D[ , indices1])[1], mean_beta1, var_beta1, var_e, type[2])
  
  # to save it more compactly?
  postvar0 = matrix(current_postvar0, length(mean_beta0), numSeq)
  postmean0 = matrix(current_postmean0, length(mean_beta0), numSeq)
  postvar1 = matrix(current_postvar1, length(mean_beta1), numSeq)
  postmean1 = matrix(current_postmean1, length(mean_beta1), numSeq)
  
  # saving it was for some evaluations anyways, not for computation of SMMED. it's fine.
  
  print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  for(t in 2:numSeq){
    
    Dt = add_MED_ms_oneatatime_data_multidm(D, y, mean_beta0, mean_beta1, 
                                             var_beta0, var_beta1, var_e,  f0, f1, type, N_seq[t], 
                                             numCandidates, k, xmin, xmax, p, alpha_seq[t], buffer_seq[t],
                                             candidates, wasserstein0, genCandidates)
    
    yt = as.vector(simulateY(Dt$addD, N_seq[t], true_beta, sigmasq, 1, true_type))
    
    # update D and y with new data
    D = c(D, Dt$addD)
    y = c(y, yt)
    
    # also save posterior means and variances
    postvar0[ , t] = diag(postvar(D, length(D), var_e, var_beta0, type[1]))
    postmean0[ , t] = postmean(y, D, length(D), mean_beta0, var_beta0, var_e, type[1])
    postvar1[ , t] = diag(postvar(D, length(D), var_e, var_beta1, type[2]))
    postmean1[ , t] = postmean(y, D, length(D), mean_beta1, var_beta1, var_e, type[2])
    
    print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
  }
  return(list("D" = D, "y" = y, "postvar0" = postvar0, "postmean0" = postmean0, 
              "postvar1" = postvar1, "postmean1" = postmean1))
}


