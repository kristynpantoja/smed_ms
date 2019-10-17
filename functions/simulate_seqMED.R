simulate_seqMED = function(true_beta, true_type, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                           f0 = NULL, f1 = NULL, type = NULL, numCandidates = 10^5, k = 4, 
                           xmin = 0, xmax = 1, p = 2, 
                           numSeq = 5, N_seq = 10, alpha_seq = NULL, buffer_seq = 0, 
                           update_prior = FALSE, wasserstein0 = 1, genCandidates = 1, seed = NULL){
  if(is.null(alpha_seq)){ # generate space-filling design for first step
    D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
                           numCandidates, k, xmin, xmax, p, alpha = 0, buffer_seq[1], 
                           genCandidates = 1, initialpt = 1)
  } else{
    D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
                           numCandidates, k, xmin, xmax, p, alpha_seq[1], buffer_seq[1], 
                           genCandidates = 1, initialpt = 1)
  }
  
  if(numSeq == 1){
    return(D1)
  }
  # otherwise, sequence:
  # some checks first:
  if(numSeq > 1 & length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
  if(numSeq > 1 & length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = c(0, ((2:numSeq) / numSeq) * (2 * p))
  if(numSeq > 1 & !is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
  if(!is.null(seed)) set.seed(seed)
  # get y1
  y1 = as.vector(simulateY(D1, N_seq[1], true_beta, sigmasq, 1, true_type))
  Nttl = sum(N_seq)
  D = D1
  y = y1
  
  current_postvar0 = diag(postvar(D, length(D), var_e, var_beta0, type[1]))
  current_postmean0 = postmean(y, D, length(D), mean_beta0, var_beta0, var_e, type[1])
  current_postvar1 = diag(postvar(D, length(D), var_e, var_beta1, type[2]))
  current_postmean1 = postmean(y, D, length(D), mean_beta1, var_beta1, var_e, type[2])
  
  postvar0 = matrix(current_postvar0, length(mean_beta0), numSeq)
  postmean0 = matrix(current_postmean0, length(mean_beta0), numSeq)
  postvar1 = matrix(current_postvar1, length(mean_beta1), numSeq)
  postmean1 = matrix(current_postmean1, length(mean_beta1), numSeq)
  
  for(t in 2:numSeq){
    if(update_prior == FALSE){
      Dt = add_MED_ms_oneatatime_data(D, y, mean_beta0, mean_beta1, 
                                      var_beta0, var_beta1, var_e,  f0, f1, type, N_seq[t], 
                                      numCandidates, k, xmin, xmax, p, alpha_seq[t], buffer_seq[t],
                                      wasserstein0, genCandidates)
    }
    if(update_prior == TRUE){
      Dt = add_MED_ms_oneatatime_data(D, y, current_postmean0, current_postmean1, 
                                      var_beta0, var_beta1, var_e,  f0, f1, type, N_seq[t], 
                                      numCandidates, k, xmin, xmax, p, alpha_seq[t], buffer_seq[t],
                                      wasserstein0, genCandidates)
    }
    
    yt = as.vector(simulateY(Dt$addD, N_seq[t], true_beta, sigmasq, 1, true_type))
    
    # update D and y with new data
    D = c(D, Dt$addD)
    y = c(y, yt)
    
    # also save posterior means and variances
    postvar0[ , t] = diag(postvar(D, length(D), var_e, var_beta0, type[1]))
    postmean0[ , t] = postmean(y, D, length(D), mean_beta0, var_beta0, var_e, type[1])
    postvar1[ , t] = diag(postvar(D, length(D), var_e, var_beta1, type[2]))
    postmean1[ , t] = postmean(y, D, length(D), mean_beta1, var_beta1, var_e, type[2])
    
    if(update_prior == TRUE){
      current_postvar0 = postvar0[ , t]
      current_postmean0 = postmean0[ , t]
      current_postvar1 = postvar1[ , t]
      current_postmean1 = postmean1[ , t]
    }
    
  }
  return(list("D" = D, "y" = y, "postvar0" = postvar0, "postmean0" = postmean0, 
              "postvar1" = postvar1, "postmean1" = postmean1))
}
