################################################################################
# simulate SeqMED, using preliminary data (or, if none, generate some)
################################################################################

SeqMED = function(
  y.in = NULL, x.in = NULL, model0, model1, error.var, 
  candidates, true.function, 
  k = 4, xmin = -1, xmax = 1, p = 1, 
  numSeq = 5, seqN = 10, alpha_seq = 1, genCandidates = 1, 
  newq = TRUE, prints = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  if(is.null(candidates)){
    if(genCandidates == 1){
      candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    } else if(genCandidates == 2){
      candidates = sort(runif(numCandidates, min = xmin, max = xmax))
    } else{
      stop("SeqMED: invalid genCandidates argument")
    }
  }
  
  # check for preliminary data
  #   if none, generate it
  if(is.null(x.in)){
    x.in = MMED(
      model0, model1, error.var, 
      seqN[1], numCandidates, k, xmin, xmax, p, alpha_seq[1],
      genCandidates = 1, initialpt = 1, var_margy0 = model0$marginalVarY, 
      var_margy1 = model1$marginalVarY)
    # x.in.idx = sample(1:length(candidates), 1)
    # x.in = candidates[x.in.idx]
  } 
  if(is.null(y.in)){
    # y.in = as.vector(simulateY(
    #   x.in, seqN[1], true_beta, error.var, 1, true_type))
    y.in = simulateY_fromfunction(x.in, true.function, error.var, 1, NULL)
  }
  x.cur = x.in
  y.cur = y.in
  x.new = c()
  y.new = c()
  
  # get posterior distributions of beta
  postbeta0 = getBetaPosterior(
    y = y.cur, X = model0$designMat(x.cur), model0$beta.mean, model0$beta.var, 
    error.var)
  postbeta1 = getBetaPosterior(
    y = y.cur, X = model1$designMat(x.cur), model1$beta.mean, model1$beta.var, 
    error.var)
  
  # postvar0.cur = postvar(x.cur, length(x.cur), error.var, beta.var0, type[1])
  # postmean0.cur = postmean(y.cur, x.cur, length(x.cur), beta.mean0, beta.var0, error.var, type[1])
  # postvar1.cur = postvar(x.cur, length(x.cur), error.var, beta.var1, type[2])
  # postmean1.cur = postmean(y.cur, x.cur, length(x.cur), beta.mean1, beta.var1, error.var, type[2])
  postvar0.cur = postbeta0$var
  postmean0.cur = postbeta0$mean
  postvar1.cur = postbeta1$var
  postmean1.cur = postbeta1$mean
  
  postvar0 = matrix(diag(postvar0.cur), length(model0$beta.mean), numSeq)
  postmean0 = matrix(postmean0.cur, length(model0$beta.mean), numSeq)
  postvar1 = matrix(diag(postvar1.cur), length(model1$beta.mean), numSeq)
  postmean1 = matrix(postmean1.cur, length(model1$beta.mean), numSeq)
  
  if(numSeq == 1){
    return(list(
      x.in = x.in, 
      y.in = y.in, 
      x.new = x.new,
      y.new = y.new,
      postvar0 = postvar0, postmean0 = postmean0, 
      postvar1 = postvar1, postmean1 = postmean1))
  }
  
  if(prints){
    print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  }
  
  # q evaluated at input points
  if(!newq){
    qs = sapply(x.cur, function(x_i) 
      q_seqmed(
        x_i, postmean0.cur, postmean1.cur, postvar0.cur, postvar1.cur, 
        error.var, type, p, alpha_seq[1]))
  }
  
  for(t in 2:numSeq){
    
    batch.idx = t - 1
    if(newq){
      Dt = SeqMED_newq_batch(
        x.cur, y.cur, model0, model1, error.var, 
        seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates,
        candidates, batch.idx)
    } else{
      Dt = SeqMED_keepq_batch(
        x.cur, y.cur, model0, model1, error.var, 
        seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates,
        candidates, batch.idx, qs = qs)
      qs = c(qs, Dt$q.new)
    }
    
    # yt = as.vector(simulateY(Dt$addD, seqN[t], true_beta, error.var, 1, true_type))
    yt = simulateY_fromfunction(
      x = Dt$addD, f = true.function, error.var = error.var, num.sims = 1, 
      seed = NULL)
    
    # update x.cur and y.cur with new data
    x.cur = c(x.cur, Dt$addD)
    y.cur = c(y.cur, yt)
    x.new = c(x.new, Dt$addD)
    y.new = c(y.new, yt)
    
    # also save posterior means and variances
    postbeta0 = getBetaPosterior(
      y = y.cur, X = model0$designMat(x.cur), model0$beta.mean, model0$beta.var, 
      error.var)
    postbeta1 = getBetaPosterior(
      y = y.cur, X = model1$designMat(x.cur), model1$beta.mean, model1$beta.var, 
      error.var)
    postvar0 = diag(postbeta0$var)
    postmean0 = postbeta0$mean
    postvar1 = diag(postbeta1$var)
    postmean1 = postbeta1$mean
    # postvar0[ , t] = diag(postvar(
    #   x.cur, length(x.cur), error.var, beta.var0, type[1]))
    # postmean0[ , t] = postmean(
    #   y.cur, x.cur, length(x.cur), beta.mean0, beta.var0, error.var, type[1])
    # postvar1[ , t] = diag(postvar(
    #   x.cur, length(x.cur), error.var, beta.var1, type[2]))
    # postmean1[ , t] = postmean(
    #   y.cur, x.cur, length(x.cur), beta.mean1, beta.var1, error.var, type[2])
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list(
    x.in = x.in, 
    y.in = y.in, 
    x.new = x.new,
    y.new = y.new,
    postvar0 = postvar0, postmean0 = postmean0, 
    postvar1 = postvar1, postmean1 = postmean1))
}




























































































# 
# SeqMEDold = function(
#   D1 = NULL, y1 = NULL, true_beta, true_type, beta.mean0, beta.mean1, 
#   beta.var0, beta.var1, error.var, f0 = NULL, f1 = NULL, type = NULL, 
#   numCandidates = 10^5, k = 4, xmin = -1, xmax = 1, p = 1, 
#   numSeq = 5, seqN = 10, alpha_seq = 1, genCandidates = 1, candidates = NULL, 
#   prints = FALSE, seed = NULL, newq = TRUE
# ){
#   if(!is.null(seed)) set.seed(seed)
#   if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
#   if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
#   if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
#   if(is.null(candidates)){
#     if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
#     if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
#   }
#   
#   # check for preliminary data
#   #   if none, generate it
#   if(is.null(D1)){
#     D1 = MMED(
#       beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type, 
#       seqN[1], numCandidates, k, xmin, xmax, p, alpha_seq[1], 
#       genCandidates = 1, initialpt = 1)
#     generatedMMED = TRUE
#   } else{
#     generatedMMED = FALSE
#   }
#   if(is.null(y1)){
#     y1 = as.vector(simulateY(D1, seqN[1], true_beta, error.var, 1, true_type))
#   }
#   Nttl = sum(seqN)
#   D = D1
#   y = y1
#   
#   current_postvar0 = postvar(D, length(D), error.var, beta.var0, type[1])
#   current_postmean0 = postmean(y, D, length(D), beta.mean0, beta.var0, error.var, type[1])
#   current_postvar1 = postvar(D, length(D), error.var, beta.var1, type[2])
#   current_postmean1 = postmean(y, D, length(D), beta.mean1, beta.var1, error.var, type[2])
#   
#   postvar0 = matrix(diag(current_postvar0), length(beta.mean0), numSeq)
#   postmean0 = matrix(current_postmean0, length(beta.mean0), numSeq)
#   postvar1 = matrix(diag(current_postvar1), length(beta.mean1), numSeq)
#   postmean1 = matrix(current_postmean1, length(beta.mean1), numSeq)
#   
#   if(numSeq == 1){
#     return(list("D" = D, "y" = y, 
#                 "postvar0" = current_postvar0, "postmean0" = current_postmean0, 
#                 "postvar1" = current_postvar1, "postmean1" = current_postmean1))
#   }
#   
#   if(prints){
#     print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
#   }
#   
#   # q evaluated at input points
#   if(!newq){
#     qs = sapply(D, function(x_i) 
#       q_seqmed(x_i, current_postmean0, current_postmean1, 
#                current_postvar0, current_postvar1, 
#                error.var, type, p, alpha_seq[1]))
#   }
#   
#   for(t in 2:numSeq){
#     
#     batch.idx = t - 1
#     if(newq){
#       Dt = SeqMED_newq_batch_old(
#         D, y, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type,
#         seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates,
#         candidates, batch.idx)
#     } else{
#       Dt = SeqMED_keepq_batch_old(
#         D, y, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type,
#         seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates,
#         candidates, batch.idx, qs = qs)
#       qs = c(qs, Dt$q.new)
#     }
#     # Dt = SeqMED_batch(
#     #   D, y, beta.mean0, beta.mean1, beta.var0, beta.var1, error.var, f0, f1, type, 
#     #   seqN[t], numCandidates, k, xmin, xmax, p, alpha_seq[t], genCandidates, 
#     #   candidates, batch.idx)
#     
#     yt = as.vector(simulateY(Dt$addD, seqN[t], true_beta, error.var, 1, true_type))
#     
#     # update D and y with new data
#     D = c(D, Dt$addD)
#     y = c(y, yt)
#     
#     # also save posterior means and variances
#     postvar0[ , t] = diag(postvar(D, length(D), error.var, beta.var0, type[1]))
#     postmean0[ , t] = postmean(y, D, length(D), beta.mean0, beta.var0, error.var, type[1])
#     postvar1[ , t] = diag(postvar(D, length(D), error.var, beta.var1, type[2]))
#     postmean1[ , t] = postmean(y, D, length(D), beta.mean1, beta.var1, error.var, type[2])
#     
#     if(prints){
#       print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
#     }
#   }
#   return(list("D" = D, "y" = y, "postvar0" = postvar0, "postmean0" = postmean0, 
#               "postvar1" = postvar1, "postmean1" = postmean1))
# }