SeqMEDvs = function(
  y.in = NULL, x.in = NULL, model0, model1, error.var, 
  candidates, true.function, true.indices, dimX, 
  k = 4, xmin = -1, xmax = 1, p = 1, 
  numSeq = 5, seqN = 10, alpha_seq = 1, 
  prints = FALSE, seed = NULL
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
  
  x.cur = x.in
  y.cur = y.in
  x.new = matrix(NA, nrow = 0, ncol = dimX)
  y.new = c()
  
  if(numSeq == 1){
    # get posterior distributions of beta
    postbeta0 = getBetaPosterior(
      y = y.cur, X = x.cur[, model0$indices, drop = FALSE], model0$beta.mean, 
      model0$beta.var, error.var)
    postbeta1 = getBetaPosterior(
      y = y.cur, X = x.cur[, model1$indices, drop = FALSE], model1$beta.mean, 
      model1$beta.var, error.var)
    
    return(list(
      x.in = x.in, 
      y.in = y.in, 
      x.new = x.new,
      y.new = y.new,
      postvar0 = diag(postbeta0$var), postmean0 = postbeta0$mean, 
      postvar1 = diag(postbeta1$var), postmean1 = postbeta1$mean))
  }
  
  for(t in 1:numSeq){
    Dt = SeqMEDvs_batch(
      initD = x.cur, inity = y.cur, model0 = model0, model1 = model1, 
      error.var = error.var, N2 = seqN[t], 
      candidates = candidates, dimX = dimX, xmin = xmin, xmax = xmax, 
      k = k, p = p, alpha = alpha_seq[t], 
      batch.idx = t)
    
    yt = simulateY_frommultivarfunction(
      x = Dt$addD, true.function = true.function, 
      error.var = error.var)
    
    # update x.cur and y.cur with new data
    x.cur = rbind(x.cur, Dt$addD)
    y.cur = c(y.cur, yt)
    x.new = rbind(x.new, Dt$addD)
    y.new = c(y.new, yt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  
  # get updated posterior distributions of beta
  #   here, x.cur is rbind(x.in, x.new) & y.cur = c(x.in, y.new)
  #   i.e. the cumulative design
  postbeta0 = getBetaPosterior(
    y = y.cur, X = x.cur[, model0$indices, drop = FALSE], model0$beta.mean, 
    model0$beta.var, error.var)
  postbeta1 = getBetaPosterior(
    y = y.cur, X = x.cur[, model1$indices, drop = FALSE], model1$beta.mean, 
    model1$beta.var, error.var)
  
  return(list(
    x.in = x.in, 
    y.in = y.in, 
    x.new = x.new,
    y.new = y.new,
    postvar0 = diag(postbeta0$var), postmean0 = postbeta0$mean, 
    postvar1 = diag(postbeta1$var), postmean1 = postbeta1$mean))
}
