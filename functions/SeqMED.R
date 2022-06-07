################################################################################
# simulate SeqMED, using preliminary data (or, if none, generate some)
################################################################################

SeqMED = function(
  y.in = NULL, x.in = NULL, model0, model1, error.var, 
  candidates, true.function, 
  xmin = -1, xmax = 1, 
  numSeq = 5, seqN = 1, alpha_seq = 1, 
  k = NULL, p = 1, genCandidates = 1, keep_trying_alpha = TRUE,
  prints = FALSE, seed = NULL, save_objectives = FALSE
){
  if(is.null(k)) k = 4 * p
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
  if(save_objectives){
    objectives_list = list()
    objectives_list[[1]] = NULL
  }
  
  # check for preliminary data
  #   if none, generate it
  if(is.null(x.in)){
    x.in = MMED(
      model0, model1, error.var, 
      seqN[1], numCandidates, k, xmin, xmax, p, alpha_seq[1],
      genCandidates = 1, initialpt = 1, var_margy0 = model0$marginalVarY, 
      var_margy1 = model1$marginalVarY)
  } 
  if(is.null(y.in)){
    y.in = simulateY_fromfunction(
      x = x.in, true.function = true.function, error.var = error.var)
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
  
  postvar0 = matrix(diag(postbeta0$var), length(model0$beta.mean), numSeq)
  postmean0 = matrix(postbeta0$mean, length(model0$beta.mean), numSeq)
  postvar1 = matrix(diag(postbeta1$var), length(model1$beta.mean), numSeq)
  postmean1 = matrix(postbeta1$mean, length(model1$beta.mean), numSeq)
  
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
  
  for(t in 2:numSeq){
    
    Dt = SeqMED_newq_batch(
      initD = x.cur, y = y.cur, model0 = model0, model1 = model1, 
      error.var = error.var, N2 = seqN[t], numCandidates = numCandidates, 
      k = k, xmin = xmin, xmax = xmax, p = p, alpha = alpha_seq[t], 
      genCandidates = genCandidates, candidates = candidates, 
      keep_trying_alpha = keep_trying_alpha, prints = prints, 
      save_objectives = save_objectives,
      batch.idx = t - 1)
    if(keep_trying_alpha){
      alpha_seq[t:numSeq] = Dt$alpha
    }
    
    yt = simulateY_fromfunction(
      x = Dt$addD, true.function = true.function, error.var = error.var)
    
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
    postvar0[, t] = diag(postbeta0$var)
    postmean0[, t] = postbeta0$mean
    postvar1[, t] = diag(postbeta1$var)
    postmean1[, t] = postbeta1$mean
    
    if(save_objectives){
      objectives_list[[t]] = Dt$objectives
    }
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  
  result = list(
    x.in = x.in, 
    y.in = y.in, 
    x.new = x.new,
    y.new = y.new,
    postvar0 = postvar0, postmean0 = postmean0, 
    postvar1 = postvar1, postmean1 = postmean1,
    alpha_seq = alpha_seq
  )
  if(save_objectives){
    result$objectives = objectives_list
  }
  return(result)
}
