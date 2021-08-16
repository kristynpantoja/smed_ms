#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

SeqMEDgp = function(
  y.in = NULL, x.in = NULL, x.in.idx = NULL, candidates, function.values, 
  xmin = 0, xmax = 1, k = 4, p = 1, 
  numSeq = 5, seqN = 3, alpha.seq = 1, buffer = 0, objective.type = 1, 
  model0 = NULL, model1 = NULL, newq = TRUE, prints = FALSE
){
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
  if(is.null(candidates)) stop("SeqMEDgp: No candidates provided!")
  
  # check preliminary data
  if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
    stop("SeqMEDgp: preliminary y.in  is given, but not corresponding x.in")
  } else if(is.null(y.in) & !is.null(x.in)){ # x is not null, y.in is null (get y.in)
    y.in = function.values[x.in.idx]
  } else if(is.null(x.in) & is.null(y.in)){ # both x.in and y.in are null, then us BH method
    stop("SeqMEDgp: need input data, at least x.in!")
  } else{
    if(length(x.in) != length(y.in)){
      stop("SeqMEDgp: length of preliminary x.in and y.in don't match!")
    }
  }
  D = x.in
  D.idx = x.in.idx
  y = y.in
  x.new = c()
  x.new.idx = c()
  y.new = c()
  
  # q evaluated at input points
  if(!newq){
    if(is.null(model0$measurement.var)){
      Kinv0 = solve(getCov(
        X1 = D, X2 = D, type = model0$type, l = model0$l, p = model0$p, 
        signal.var = model0$signal.var))
    } else{
      Kinv0 = solve(getCov(
        X1 = D, X2 = D, type = model0$type, l = model0$l, p = model0$p, 
        signal.var = model0$signal.var) + 
          model0$measurement.var * diag(length(D)))
    }
    if(is.null(model1$measurement.var)){
      Kinv1 = solve(getCov(
        X1 = D, X2 = D, type = model1$type, l = model1$l, p = model1$p, 
        signal.var = model1$signal.var))
    } else{
      Kinv1 = solve(getCov(
        X1 = D, X2 = D, type = model1$type, l = model1$l, p = model1$p, 
        signal.var = model1$signal.var) + 
          model1$measurement.var * diag(length(D)))
    }
    qs = rep(NA, length(D))
    if(!(objective.type %in% c(0, 1, 3, 4, 5))){
      stop("SeqMEDgp: to keep q, need objective.type == 1, 3, 4, or 5")
    } else{
      if(objective.type == 1){ # buffer
        qs = sapply(D, function(x_i) 
          q_gp(
            x_i, Kinv0, Kinv1, D, y, p, alpha.seq[1], buffer, model0, model1))
      }
      if(objective.type %in% c(0, 3, 5)){
        qs = rep(1, length(D))
      }
      if(objective.type == 4){ # cap q
        qs = sapply(D, function(x_i) 
          qcap_gp(x_i, Kinv0, Kinv1, D, y, p, alpha.seq[1], model0, model1))
      }
    }
  }
  
  for(t in 1:numSeq){
    
    batch.idx = t
    
    if(newq){
      Dt = SeqMEDgp_newq_batch(
        initD = D, y = y, N2 = seqN[t], numCandidates = numCandidates, k = k, p = p, 
        xmin = xmin, xmax = xmax, alpha = alpha.seq[t], candidates = candidates, 
        batch.idx = batch.idx, buffer = buffer, objective.type = objective.type, 
        model0 = model0, model1 = model1)
    } else{
      Dt = SeqMEDgp_keepq_batch(
        initD = D, y = y, N2 = seqN[t], numCandidates = numCandidates, k = k, p = p, 
        xmin = xmin, xmax = xmax, alpha = alpha.seq[t], candidates = candidates, 
        batch.idx = batch.idx, buffer = buffer, objective.type = objective.type,
        model0 = model0, model1 = model1, qs = qs)
      qs = c(qs, Dt$q.new)
    }
    
    yt = function.values[Dt$indices]
    
    # update D and y with new data
    D = c(D, Dt$addD)
    D.idx = c(D.idx, Dt$indices)
    y = c(y, yt)
    x.new = c(x.new, Dt$addD)
    x.new.idx = c(x.new.idx, Dt$indices)
    y.new = c(y.new, yt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list(
    x = x.in, 
    x.idx = x.in.idx, 
    y = y.in, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    function.values = function.values
  ))
}
