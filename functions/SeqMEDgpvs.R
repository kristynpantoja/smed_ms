########################################################################
### SeqMED GP for variable selection, one-at-a-time greedy algorithm ###
########################################################################

SeqMEDgpvs = function(
  y.in = NULL, x.in = NULL, x.in.idx = NULL,
  candidates, function.values, 
  xmin = 0, xmax = 1, k = 4, p = 1, 
  numSeq = 5, seqN = 3, alpha.seq = 1, buffer = 0, objective.type = 1, 
  model0 = NULL, model1 = NULL, newq = TRUE, prints = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
  if(is.null(candidates)) stop("SeqMEDgpvs: No candidates provided!")
  
  # check preliminary data
  if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
    stop("SeqMEDgpvs: preliminary y.in  is given, but not corresponding x.in")
  } else if(is.null(y.in) & !is.null(x.in)){ # x is not null, y.in is null (get y.in)
    y.in = function.values[x.in.idx]
  } else if(is.null(x.in) & is.null(y.in)){ # both x.in and y.in are null, then us BH method
    stop("SeqMEDgpvs: need input data, at least x.in!")
  } else{
    if(nrow(x.in) != length(y.in)){
      stop("SeqMEDgpvs: length of preliminary x.in and y.in don't match!")
    }
  }
  x.cur = x.in
  y.cur = y.in
  x.new = c()
  x.new.idx = c()
  y.new = c()
  
  # q evaluated at input points
  if(!newq){
    initD0 = x.cur[ , model0$indices, drop = FALSE]
    initD1 = x.cur[ , model1$indices, drop = FALSE]
    if(is.null(model0$measurement.var)){
      Kinv0 = solve(getCov(
        X1 = initD0, X2 = initD0, 
        type = model0$type, l = model0$l, p = model0$p, 
        signal.var = model0$signal.var))
    } else{
      Kinv0 = solve(getCov(
        X1 = initD0, X2 = initD0, 
        type = model0$type, l = model0$l, p = model0$p, 
        signal.var = model0$signal.var) + 
          model0$measurement.var * diag(nrow(x.cur)))
    }
    if(is.null(model1$measurement.var)){
      Kinv1 = solve(getCov(
        X1 = initD1, X2 = initD1, 
        type = model1$type, l = model1$l, p = model1$p, 
        signal.var = model1$signal.var))
    } else{
      Kinv1 = solve(getCov(
        X1 = initD1, X2 = initD1, 
        type = model1$type, l = model1$l, p = model1$p, 
        signal.var = model1$signal.var) + 
          model1$measurement.var * diag(nrow(x.cur)))
    }
    qs = rep(NA, nrow(x.cur))
    if(!(objective.type %in% c(0, 1, 3, 4, 5))){
      stop("SeqMEDgpvs: to keep q, need objective.type == 1, 3, 4, or 5")
    } else{
      if(objective.type == 1){ # buffer
        qs = apply(x.cur, 1, function(x_i) 
          q_gpvs(
            x_i, Kinv0, Kinv1, initD0, initD1, y.cur, p, alpha.seq[1], buffer, 
            model0, model1))
      }
      if(objective.type %in% c(0, 3, 5)){
        qs = rep(1, nrow(x.cur))
      }
      if(objective.type == 4){ # cap q
        qs = apply(x.cur, 1, function(x_i) 
          qcap_gpvs(
            x_i, Kinv0, Kinv1, initD0, initD1, y.cur, p, alpha.seq[1], 
            model0, model1))
      }
    }
  }
  
  for(t in 1:numSeq){
    
    batch.idx = t
    
    if(newq){
      Dt = SeqMEDgpvs_newq_batch(
        initD = x.cur, y = y.cur, N2 = seqN[t], numCandidates = numCandidates, 
        k = k, p = p,
        xmin = xmin, xmax = xmax, alpha = alpha.seq[t], candidates = candidates, 
        batch.idx = batch.idx, buffer = buffer, objective.type = objective.type, 
        model0 = model0, model1 = model1)
      
    } else{
      Dt = SeqMEDgpvs_keepq_batch(
        initD = x.cur, y = y.cur, N2 = seqN[t], numCandidates = numCandidates, 
        k = k, p = p,
        xmin = xmin, xmax = xmax, alpha = alpha.seq[t], candidates = candidates, 
        batch.idx = batch.idx, buffer = buffer, objective.type = objective.type, 
        model0 = model0, model1 = model1, qs = qs)
      qs = c(qs, Dt$q.new)
    }
    
    yt = function.values[Dt$indices]
    
    # update x.cur and y.cur with new data
    x.cur = rbind(x.cur, Dt$addD)
    y.cur = c(y.cur, yt)
    x.new = rbind(x.new, Dt$addD)
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
