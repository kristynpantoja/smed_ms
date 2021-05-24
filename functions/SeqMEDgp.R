#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

SeqMEDgp = function(
  y0 = NULL, x0 = NULL, x0.idx = NULL, candidates, function.values, 
  xmin = 0, xmax = 1, k = 4, p = 1, 
  numSeq = 5, seqN = 3, alpha.seq = 1, buffer = 0, objective.type = 1, 
  init.as.stage = FALSE, prints = FALSE, seed = NULL, 
  model0 = NULL, model1 = NULL, noise = TRUE, measurement.var = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
  if(is.null(candidates)) stop("SeqMEDgp: No candidates provided!")
  
  # check preliminary data
  if(is.null(x0) & !is.null(y0)){ # x0 is null, y0 is not null
    stop("SeqMEDgp : preliminary y0  is given, but not corresponding x0")
  } else if(is.null(y0) & !is.null(x0)){ # x is not null, y0 is null (get y0)
    y0 = function.values[x0.idx]
  } else if(is.null(x0) & is.null(y0)){ # both x0 and y0 are null, then us BH method
    stop("SeqMEDgp: need input data, at least x0!")
  } else{
    if(length(x0) != length(y0)){
      stop("SeqMEDgp : length of preliminary x0 and y0 don't match!")
    }
  }
  D = x0
  D.idx = x0.idx
  y = y0
  x.new = c()
  x.new.idx = c()
  y.new = c()
  
  if(init.as.stage){ # if initial data is its own stage
    if(numSeq == 1){
      return(list(
        x = x0, 
        x.idx = x0.idx, 
        y = y0, 
        x.new = x.new,
        x.new.idx = x.new.idx,
        y.new = y.new,
        function.values = function.values, 
        # old outputs, in case they're needed
        D = D, 
        D.idx = D.idx, 
        y = y
      ))
    }
    
    if(prints){
      print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
    }
    tStart = 2
  } else{
    tStart = 1
  }
  
  for(t in tStart:numSeq){
    
    if(tStart == 2){
      batch.idx = t - 1
    } else{
      batch.idx = t
    }
    Dt = SeqMEDgp_batch(
      initD = D, y = y, N2 = seqN[t],
      k = k, p = p, xmin = xmin, xmax = xmax,
      alpha = alpha.seq[t], candidates = candidates, batch.idx = batch.idx, 
      buffer = buffer, objective.type = objective.type, 
      model0 = model0, model1 = model1)
    
    yt = function.values[Dt$indices]
    if(noise){
      if(is.null(measurement.var)){
        stop("SeqMEDgp: noise = TRUE, but measurement.var = NULL.")
      } else{
        yt = yt + rnorm(seqN[t], 0, sqrt(measurement.var))
      }
    }
    
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
    x = x0, 
    x.idx = x0.idx, 
    y = y0, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    function.values = function.values, 
    # old outputs, in case they're needed
    D = D, 
    D.idx = D.idx, 
    y = y
  ))
}

##################################################
### SeqMED GP for Variable Selection, 1D vs 2D ###
##################################################

# variable selection case
# generate_SMMEDgpvs
SeqMEDgpvs = function(
  # true_y, type_true = NULL, l_true = NULL, indices_true = NULL, 
  y0, x0 = NULL, x0.idx = NULL, function.values, 
  type = c(1, 1), l = c(0.1, 0.1), idx0, idx1, signal.var = 1, N2 = 11, 
  candidates, k = 4, p = 1, 
  xmin = 0, xmax = 1, nugget = NULL, 
  numSeq = 5, seqN = 3, alpha.seq = 1, buffer_seq = 0, 
  numDims = NULL, 
  seed = NULL, algorithm = 1, prints = FALSE, init.as.stage = TRUE
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
  if(is.null(candidates)) stop("SeqMEDgp: No candidates provided!")
  
  # # check candidates
  # if(is.null(candidates)){
  #   if(!is.null(numDims)){
  #     xseq = seq(from = xmin, to = xmax, length.out = numCandidates)
  #     xseqlist = list()
  #     for(i in 1:numDims) xseqlist[[i]] = xseq 
  #     candidates = as.matrix(expand.grid(xseqlist))
  #   } else{
  #     stop("SeqMEDgpvs: Candidates and numDims == NULL. Cannot generate 
  #     canddiates without numDims.")
  #   }
  # }
  # check preliminary data
  if(is.null(x0) & !is.null(y0)){ # x0 is null, y0 is not null
    stop("SeqMEDgp : preliminary y0  is given, but not corresponding x0")
  } else if(is.null(y0) & !is.null(x0)){ # x is not null, y0 is null (get y0)
    y0 = function.values[x0.idx]
  } else if(is.null(x0) & is.null(y0)){ # both x0 and y0 are null, then us BH method
    stop("SeqMEDgp: need input data, at least x0!")
  } else{
    if(length(x0) != length(y0)){
      stop("SeqMEDgp : length of preliminary x0 and y0 don't match!")
    }
  }
  Nttl = sum(seqN)
  D = x0
  D.idx = x0.idx
  y = y0
  x.new = c()
  x.new.idx = c()
  y.new = c()
  
  # # check initD and inity
  # if(is.null(x0)){
  #   if(is.null(initN)){
  #     warning("initN and initD not given; going to generate seqN[1] initial points")
  #     initN = seqN[1]
  #   }
  #   x0.idx = sample(1:dim(candidates)[1], initN)
  #   x0 = candidates[initD_indices, ]
  # }
  # # inity
  # if(is.null(inity)){
  #   inity = true_y[initD_indices]
  # } else{
  #   y1 = inity
  # }
  # D1 = x0
  
  if(init.as.stage){ # if initial data is its own stage
    if(numSeq == 1){
      return(list(
        x = x0, 
        x.idx = x0.idx, 
        y = y0, 
        x.new = x.new,
        x.new.idx = x.new.idx,
        y.new = y.new,
        function.values = function.values, 
        # old outputs, in case they're needed
        D = D, 
        D.idx = D.idx, 
        y = y
      ))
    }
    
    if(prints){
      print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
    }
    tStart = 2
  } else{
    tStart = 1
  }
  
  for(t in tStart:numSeq){
    
    if(tStart == 2){
      batch.idx = t - 1
    } else{
      batch.idx = t
    }
    Dt = SeqMEDgp_batch(
      initD = D, y = y, type = type, l = l, signal.var = signal.var, N2 = seqN[t],
      k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget, 
      alpha = alpha.seq[t], candidates = candidates, batch.idx = batch.idx, 
      obj_fn = obj_fn)
    
    yt = function.values[Dt$indices]
    
    
    if(algorithm == 1){
      Dt = add_MMEDgpvs(initD = D, y = y, type = type_hypotheses, l = l_hypotheses, 
                        indices0 = indices0, indices1 = indices1, var_e = var_e, 
                        N2 = seqN[t], numCandidates = NULL, k = k, p = p, 
                        xmin = xmin, xmax = xmax, nugget = nugget, 
                        alpha = alpha.seq[t], buffer = buffer_seq[t], 
                        genCandidates = genCandidates, candidates = candidates, 
                        algorithm = algorithm, batch.idx = batch.idx)
    }
    yt = function.values[Dt$indices]
    
    # update D and y with new data
    D = rbind(D, Dt$addD)
    D.idx = c(D.idx, Dt$indices)
    y = c(y, yt)
    x.new = rbind(x.new, Dt$addD)
    x.new.idx = c(x.new.idx, Dt$indices)
    y.new = c(y.new, yt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list(
    x = x0, 
    x.idx = x0.idx, 
    y = y0, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    function.values = function.values, 
    # old outputs, in case they're needed
    D = D, 
    D.idx = D.idx, 
    y = y
  ))
  
  
  return(list(
    "initD" = initD, 
              "addD" = D, 
              "D" = c(initD, D),
              "candidates" = candidates, 
              "indices" = indices
    ))
}
