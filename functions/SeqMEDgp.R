#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

SeqMEDgp = function(
  y.in = NULL, x.in = NULL, x.in.idx = NULL, candidates, function.values, 
  xmin = -1, xmax = 1, k = NULL, p = 1, 
  n = NULL, numSeq = 5, seqN = 1, alpha.seq = 1, 
  model0 = NULL, model1 = NULL, prints = FALSE
){
  if(is.null(k)) k = 4 * p
  if(!is.null(n)){
    numSeq = n
    seqN = 1
  }
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
  if(is.null(candidates)) stop("SeqMEDgp: No candidates provided!")
  
  # check preliminary data
  if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
    stop("SeqMEDgp: preliminary y.in  is given, but not corresponding x.in")
  } else if(is.null(x.in)){
    x.in.idx = ceiling(numx / 2)
    x.in = candidates[x.in.idx]
  }
  if(is.null(y.in)){ # x is not null, y.in is null (get y.in)
    y.in = function.values[x.in.idx]
  } else{
    if(length(x.in) != length(y.in)){
      stop("SeqMEDgp: length of preliminary x.in and y.in don't match!")
    }
  }
  x.cur = x.in
  y.cur = y.in
  x.new = c()
  x.new.idx = c()
  y.new = c()
  
  if(numSeq == 1){
    return(list(
      x.in = x.in, 
      x.in.idx = x.in.idx, 
      y.in = y.in, 
      x.new = x.new,
      x.new.idx = x.new.idx, 
      y.new = y.new,
      function.values = function.values))
  }
  
  if(prints){
    print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  }
  
  for(t in 2:numSeq){
    batch.idx = t - 1
    
    Dt = SeqMEDgp_newq_batch(
      initD = x.cur, y = y.cur, N2 = seqN[t], numCandidates = numCandidates, 
      k = k, p = p, 
      xmin = xmin, xmax = xmax, alpha = alpha.seq[t], candidates = candidates, 
      batch.idx = batch.idx, model0 = model0, model1 = model1)
    
    yt = function.values[Dt$indices]
    
    # update x.cur and y.cur with new data
    x.cur = c(x.cur, Dt$addD)
    y.cur = c(y.cur, yt)
    x.new = c(x.new, Dt$addD)
    x.new.idx = c(x.new.idx, Dt$indices)
    y.new = c(y.new, yt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list(
    x.in = x.in, 
    x.in.idx = x.in.idx, 
    y.in = y.in, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    function.values = function.values
  ))
}
