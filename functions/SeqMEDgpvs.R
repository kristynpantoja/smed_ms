########################################################################
### SeqMED GP for variable selection, one-at-a-time greedy algorithm ###
########################################################################

# NOTES:
# new arguments (compared to SeqMEDgp):
#   seed = NULL
# old arguments that need to be re-defined: p = # dimensions

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
            x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha.seq[1], buffer, 
            model0, model1))
      }
      if(objective.type %in% c(0, 3, 5)){
        qs = rep(1, nrow(x.cur))
      }
      if(objective.type == 4){ # cap q
        qs = apply(x.cur, 1, function(x_i) 
          qcap_gpvs(
            x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha.seq[1], 
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

# ##################################################
# ### SeqMED GP for Variable Selection, 1D vs 2D ###
# ##################################################
# 
# # variable selection case
# # generate_SMMEDgpvs
# SeqMEDgpvs = function(
#   # true_y, type_true = NULL, l_true = NULL, indices_true = NULL, 
#   y.in, x.in = NULL, x.in.idx = NULL, function.values, 
#   type = c(1, 1), l = c(0.1, 0.1), idx.in, idx1, signal.var = 1, N2 = 11, 
#   candidates, k = 4, p = 1, 
#   xmin = 0, xmax = 1, nugget = NULL, 
#   numSeq = 5, seqN = 3, alpha.seq = 1, buffer_seq = 0, 
#   numDims = NULL, 
#   seed = NULL, algorithm = 1, prints = FALSE, init.as.stage = TRUE
# ){
#   if(!is.null(seed)) set.seed(seed)
#   if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
#   if(numSeq > 1 & is.null(alpha.seq)) alpha.seq = rep(1, numSeq)
#   if(numSeq > 1 & length(alpha.seq) == 1) alpha.seq = rep(alpha.seq, numSeq)
#   if(is.null(candidates)) stop("SeqMEDgp: No candidates provided!")
#   
#   # # check candidates
#   # if(is.null(candidates)){
#   #   if(!is.null(numDims)){
#   #     xseq = seq(from = xmin, to = xmax, length.out = numCandidates)
#   #     xseqlist = list()
#   #     for(i in 1:numDims) xseqlist[[i]] = xseq 
#   #     candidates = as.matrix(expand.grid(xseqlist))
#   #   } else{
#   #     stop("SeqMEDgpvs: Candidates and numDims == NULL. Cannot generate 
#   #     canddiates without numDims.")
#   #   }
#   # }
#   # check preliminary data
#   if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
#     stop("SeqMEDgp : preliminary y.in  is given, but not corresponding x.in")
#   } else if(is.null(y.in) & !is.null(x.in)){ # x is not null, y.in is null (get y.in)
#     y.in = function.values[x.in.idx]
#   } else if(is.null(x.in) & is.null(y.in)){ # both x.in and y.in are null, then us BH method
#     stop("SeqMEDgp: need input data, at least x.in!")
#   } else{
#     if(length(x.in) != length(y.in)){
#       stop("SeqMEDgp : length of preliminary x.in and y.in don't match!")
#     }
#   }
#   Nttl = sum(seqN)
#   D = x.in
#   D.idx = x.in.idx
#   y = y.in
#   x.new = c()
#   x.new.idx = c()
#   y.new = c()
#   
#   # # check initD and inity
#   # if(is.null(x.in)){
#   #   if(is.null(initN)){
#   #     warning("initN and initD not given; going to generate seqN[1] initial points")
#   #     initN = seqN[1]
#   #   }
#   #   x.in.idx = sample(1:dim(candidates)[1], initN)
#   #   x.in = candidates[initD_indices, ]
#   # }
#   # # inity
#   # if(is.null(inity)){
#   #   inity = true_y[initD_indices]
#   # } else{
#   #   y1 = inity
#   # }
#   # D1 = x.in
#   
#   if(init.as.stage){ # if initial data is its own stage
#     if(numSeq == 1){
#       return(list(
#         x = x.in, 
#         x.idx = x.in.idx, 
#         y = y.in, 
#         x.new = x.new,
#         x.new.idx = x.new.idx,
#         y.new = y.new,
#         function.values = function.values, 
#         # old outputs, in case they're needed
#         D = D, 
#         D.idx = D.idx, 
#         y = y
#       ))
#     }
#     
#     if(prints){
#       print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
#     }
#     tStart = 2
#   } else{
#     tStart = 1
#   }
#   
#   for(t in tStart:numSeq){
#     
#     if(tStart == 2){
#       batch.idx = t - 1
#     } else{
#       batch.idx = t
#     }
#     Dt = SeqMEDgp_batch(
#       initD = D, y = y, type = type, l = l, signal.var = signal.var, N2 = seqN[t],
#       k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget, 
#       alpha = alpha.seq[t], candidates = candidates, batch.idx = batch.idx, 
#       obj_fn = obj_fn)
#     
#     yt = function.values[Dt$indices]
#     
#     
#     if(algorithm == 1){
#       Dt = add_MMEDgpvs(initD = D, y = y, type = type_hypotheses, l = l_hypotheses, 
#                         indices0 = indices0, indices1 = indices1, var_e = var_e, 
#                         N2 = seqN[t], numCandidates = NULL, k = k, p = p, 
#                         xmin = xmin, xmax = xmax, nugget = nugget, 
#                         alpha = alpha.seq[t], buffer = buffer_seq[t], 
#                         genCandidates = genCandidates, candidates = candidates, 
#                         algorithm = algorithm, batch.idx = batch.idx)
#     }
#     yt = function.values[Dt$indices]
#     
#     # update D and y with new data
#     D = rbind(D, Dt$addD)
#     D.idx = c(D.idx, Dt$indices)
#     y = c(y, yt)
#     x.new = rbind(x.new, Dt$addD)
#     x.new.idx = c(x.new.idx, Dt$indices)
#     y.new = c(y.new, yt)
#     
#     if(prints){
#       print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
#     }
#   }
#   return(list(
#     x = x.in, 
#     x.idx = x.in.idx, 
#     y = y.in, 
#     x.new = x.new,
#     x.new.idx = x.new.idx,
#     y.new = y.new,
#     function.values = function.values, 
#     # old outputs, in case they're needed
#     D = D, 
#     D.idx = D.idx, 
#     y = y
#   ))
#   
#   
#   return(list(
#     "initD" = initD, 
#               "addD" = D, 
#               "D" = c(initD, D),
#               "candidates" = candidates, 
#               "indices" = indices
#     ))
# }
