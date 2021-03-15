SeqMEDgp = function(
  y0 = NULL, x0 = NULL, x0.idx = NULL, candidates, function.values, nugget = NULL,
  type, l, error.var = 1, xmin = 0, xmax = 1, k = 4, p = 1, numSeq = 5, seqN = 3, 
  alpha_seq = 1, prints = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
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
  
  if(numSeq == 1){
    return(list("D" = D, "D.idx" = D.idx, "y" = y))
  }
  
  if(prints){
    print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  }
  for(t in 2:numSeq){
    
    batch.idx = t - 1
    Dt = SeqMEDgp_batch(
      initD = D, y = y, type = type, l = l, error.var = error.var, N2 = seqN[t],
      k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget, 
      alpha = alpha_seq[t], candidates = candidates, batch.idx = batch.idx)
    
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



# variable selection case
# generate_SMMEDgpvs
SeqMEDgpvs = function(
  true_y, type_true = NULL, l_true = NULL, 
  indices_true = NULL, initD = NULL, initD_indices = NULL, 
  inity = NULL, 
  type_hypotheses = c(1, 1), l_hypotheses = c(0.1, 0.1), 
  indices0, indices1, var_e = 1, N2 = 11, 
  numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, nugget = NULL, 
  numSeq = 5, N_seq = 3, alpha_seq = 1, buffer_seq = 0, 
  genCandidates = 1, candidates = NULL, numDims = NULL, 
  seed = NULL, algorithm = 1, prints = FALSE
){
  # (type_true = NULL, l_true = NULL, 
  #  indices_true = NULL, true_y = NULL, 
  #  type_hypotheses, l_hypotheses, indices0, indices1, 
  #  var_e = 1, xmin = 0, xmax = 1, nugget = 1e-1, numCandidates = NULL, k = 4, p = 1, 
  #  initN = NULL, initD = NULL, initD_indices = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
  #  alpha_seq = 1, buffer_seq = 0, numDimsMax = NULL, 
  #  genCandidates = 1, candidates = NULL, algorithm = 1, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  
  # check candidates
  if(is.null(candidates)){
    if(!is.null(numDims)){
      xseq = seq(from = xmin, to = xmax, length.out = numCandidates)
      xseqlist = list()
      for(i in 1:numDims) xseqlist[[i]] = xseq 
      candidates = as.matrix(expand.grid(xseqlist))
    }
  }
  
  # check sequential settings
  if(length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
  if(length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(!is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
  # check initD and inity
  if(is.null(initD)){
    if(is.null(initN)){
      warning("initN and initD not given; going to generate N_seq[1] initial points")
      initN = N_seq[1]
    }
    initD_indices = sample(1:dim(candidates)[1], initN)
    initD = candidates[initD_indices, ]
  }
  # inity
  if(is.null(inity)){
    inity = true_y[initD_indices]
  } else{
    y1 = inity
  }
  D1 = initD
  
  if(numSeq == 1){
    return(list("D" = D1,
                "indices" = initD_indices,
                "y" = y1,
                "candidates" = candidates))
  }
  y = y1
  D = D1
  indices = c()
  if(prints){
    print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  }
  
  for(t in 1:numSeq){
    
    batch.idx = t - 1
    Dt = add_MMEDgpvs(initD = D, y = y, type = type_hypotheses, l = l_hypotheses, 
                      indices0 = indices0, indices1 = indices1, var_e = var_e, 
                      N2 = N_seq[t], numCandidates = NULL, k = k, p = p, 
                      xmin = xmin, xmax = xmax, nugget = nugget, 
                      alpha = alpha_seq[t], buffer = buffer_seq[t], 
                      genCandidates = genCandidates, candidates = candidates, 
                      algorithm = algorithm, batch.idx = batch.idx)
    
    indt = Dt$indices
    yt = true_y[indt]
    
    # update D and y with new data
    D = rbind(D, Dt$addD)
    y = c(y, yt)
    indices = c(indices, indt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  return(list(
    "initD" = initD, 
              "addD" = D, 
              "D" = c(initD, D),
              "candidates" = candidates, 
              "indices" = indices
    ))
}
