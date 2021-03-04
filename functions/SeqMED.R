################################################################################
# simulate SeqMED, using preliminary data (or, if none, generate some)
################################################################################

# simulate_seqMED
SeqMED = function(
  D1 = NULL, y1 = NULL, true_beta, true_type, beta.mean0, beta.mean1, 
  beta.var0, beta.var1, error.var, f0 = NULL, f1 = NULL, type = NULL, 
  numCandidates = 10^5, k = 4, xmin = -1, xmax = 1, p = 1, 
  numSeq = 5, seqN = 10, alpha_seq = 1, genCandidates = 1, candidates = NULL, 
  prints = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  # some checks
  if(is.null(D1)){
    D1 = MMED(
      beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type, 
      seqN[1], numCandidates, k, xmin, xmax, p, alpha_seq[1], 
      genCandidates = 1, initialpt = 1)
    generatedMMED = TRUE
  } else{
    generatedMMED = FALSE
  }
  
  # get y1
  if(is.null(y1)){
    y1 = as.vector(simulateY(D1, seqN[1], true_beta, error.var, 1, true_type))
  }
  Nttl = sum(seqN)
  D = D1
  y = y1
  
  current_postvar0 = diag(postvar(D, length(D), error.var, beta.var0, type[1]))
  current_postmean0 = postmean(y, D, length(D), beta.mean0, beta.var0, error.var, type[1])
  current_postvar1 = diag(postvar(D, length(D), error.var, beta.var1, type[2]))
  current_postmean1 = postmean(y, D, length(D), beta.mean1, beta.var1, error.var, type[2])
  
  postvar0 = matrix(current_postvar0, length(beta.mean0), numSeq)
  postmean0 = matrix(current_postmean0, length(beta.mean0), numSeq)
  postvar1 = matrix(current_postvar1, length(beta.mean1), numSeq)
  postmean1 = matrix(current_postmean1, length(beta.mean1), numSeq)
  
  if(numSeq == 1){
    return(list("D" = D, "y" = y, 
                "postvar0" = current_postvar0, "postmean0" = current_postmean0, 
                "postvar1" = current_postvar1, "postmean1" = current_postmean1))
  }
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  if(prints){
    print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  }
  for(t in 2:numSeq){
    
    batch.idx = t - 1
    Dt = SeqMED_batch(
      D, y, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type, 
      seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates, 
      candidates, batch.idx)
    
    yt = as.vector(simulateY(Dt$addD, seqN[t], true_beta, error.var, 1, true_type))
    
    # update D and y with new data
    D = c(D, Dt$addD)
    y = c(y, yt)
    
    # also save posterior means and variances
    postvar0[ , t] = diag(postvar(D, length(D), error.var, beta.var0, type[1]))
    postmean0[ , t] = postmean(y, D, length(D), beta.mean0, beta.var0, error.var, type[1])
    postvar1[ , t] = diag(postvar(D, length(D), error.var, beta.var1, type[2]))
    postmean1[ , t] = postmean(y, D, length(D), beta.mean1, beta.var1, error.var, type[2])
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list("D" = D, "y" = y, "postvar0" = postvar0, "postmean0" = postmean0, 
              "postvar1" = postvar1, "postmean1" = postmean1))
}