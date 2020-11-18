# for multidimensional case... like LASSO or ridge
generate_SMMEDvs = function(mean_beta_full, beta_true = NULL, indices_true = NULL, 
                            indices0, indices1, mean_beta0 = NULL, mean_beta1 = NULL, 
                            var_e = 1, var_beta = NULL, var_beta0 = NULL, var_beta1 = NULL, 
                            xmin = -1, xmax = 1, numCandidates = 10^5, k = 4, p = 1,
                            initD = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
                            alpha_seq = NULL, buffer_seq = 0, candidates = NULL, 
                            wasserstein0 = 1, genCandidates = 1, seed = NULL, 
                            algorithm = 1){
  
  # first check if some things are null
  # some checks!!!!!!!!!!!!!!!!!!!!!!
  if(is.null(mean_beta_full)) stop("mean_beta_full is not given")
  if(is.null(beta_true)){
    if(is.null(indices_true)){
      stop("neither beta_true nor indices_true are given")
    } else{
      beta_true = mean_beta_full[indices_true]
    }
  } else{
    if(!is.null(indices_true)){
      if(!all.equal(beta_true, mean_beta_full[indices_true])){
        stop("beta_true does not match indices_true")
      }
    }
  }
  p0 = length(mean_beta0)
  p1 = length(mean_beta1)
  if(is.null(var_beta0) & is.null(var_beta1)){
    if(is.null(var_beta)){
      stop("neither var_beta0, var_beta1, nor var_beta are given")
    } else{
      var_beta0 = diag(rep(var_beta, p0))
      var_beta1 = diag(rep(var_beta, p1))
    }
  }
  if(is.null(initD)){
    stop("no initial data! haven't written this case for multiple dimensions yet. going to just randomly select points.")
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
  # sequence checks
  if(length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
  if(length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(!is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
  if(!is.null(seed)) set.seed(seed)
  # get y1 (if not given)
  if(is.null(inity)){
    y1 = as.vector(simulateYvs(D1[ , indices_true], N_seq[1], beta_true, sigmasq, 1, seed = seed))
  } else{
    y1 = inity
  }
  y = y1
  D = D1
  if(numSeq == 1){
    return(list("D" = D1,
                "y" = y1))
  }
  
  # why do I do this? why only take the diagonal of postvar?
  current_postvar0 = diag(postvar(D[ , indices0], initN, var_e, var_beta0, type = NULL))
  current_postmean0 = postmean(y, D[ , indices0], initN, mean_beta0, var_beta0, var_e, type = NULL)
  current_postvar1 = diag(postvar(D[ , indices1], initN, var_e, var_beta1, type = NULL))
  current_postmean1 = postmean(y, D[ , indices1], initN, mean_beta1, var_beta1, var_e, type = NULL)
  
  # to save it more compactly; only need pointwise variance, no need to save whole covariance matrix
  # first one is for initial data
  # the other numSeq - 1 ones will be replaced during the numSeq loop
  postvar0 = matrix(current_postvar0, length(mean_beta0), numSeq)
  postmean0 = matrix(current_postmean0, length(mean_beta0), numSeq)
  postvar1 = matrix(current_postvar1, length(mean_beta1), numSeq)
  postmean1 = matrix(current_postmean1, length(mean_beta1), numSeq)
  
  # get candidates
  pfull = length(mean_beta_full)
  if(is.null(candidates)){
    candidates = matrix(runif(n = pfull * numCandidates, min = xmin, max = xmax), 
                        nrow = numCandidates, ncol = pfull)
  }
  
  print(paste("finished ", 0, " out of ", numSeq, " steps", sep = ""))
  for(t in 1:numSeq){
    seed = seed + t
    
    Dt = add_MMEDvs(initD = D, inity = y, mean_beta_full, beta_true, indices_true, 
                    indices0, indices1, mean_beta0, mean_beta1, 
                    var_e, var_beta, var_beta0, var_beta1,
                    N2 = N_seq[t], xmin, xmax, numCandidates, k, p, 
                    alpha = alpha_seq[t], buffer = buffer_seq[t], candidates,
                    wasserstein0, genCandidate, algorithm = algorithm)
    
    yt = as.vector(simulateYvs(Dt$addD[ , indices_true], N_seq[t], beta_true, sigmasq, 1, seed = seed))
    
    # update D and y with new data
    D = rbind(D, Dt$addD)
    y = c(y, yt)
    
    # also save posterior means and variances
    currentN = dim(D)[1]
    postvar0[ , t] = diag(postvar(D[ , indices0], currentN, var_e, var_beta0, type = NULL))
    postmean0[ , t] = postmean(y, D[ , indices0], currentN, mean_beta0, var_beta0, var_e, type = NULL)
    postvar1[ , t] = diag(postvar(D[ , indices1], currentN, var_e, var_beta1, type = NULL))
    postmean1[ , t] = postmean(y, D[ , indices1], currentN, mean_beta1, var_beta1, var_e, type = NULL)
    
    print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
  }
  return(list("initD" = initD, "addD" = D[(initN + 1):dim(D)[1], ], "D" = D, "y" = y, 
              "postvar0" = postvar0, "postmean0" = postmean0, 
              "postvar1" = postvar1, "postmean1" = postmean1))
}


