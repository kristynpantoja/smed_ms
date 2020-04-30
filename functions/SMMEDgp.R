# variable selection case

generate_SMMEDgpvs = function(true_y, type_true = NULL, l_true = NULL, 
                              indices_true = NULL, initD = NULL, initD_indices = NULL, 
                              inity = NULL, 
                              type_hypotheses = c(1, 1), l_hypotheses = c(0.1, 0.1), 
                              indices0, indices1, var_e = 1, N2 = 11, 
                              numCandidates = 10^5, k = 4, p = 1, 
                              xmin = 0, xmax = 1, nugget = NULL, 
                              numSeq = 5, N_seq = 3, alpha_seq = 1, buffer_seq = 0, 
                              genCandidates = 1, candidates = NULL, numDims = NULL, 
                              seed = NULL, algorithm = 1){
# (type_true = NULL, l_true = NULL, 
#  indices_true = NULL, true_y = NULL, 
#  type_hypotheses, l_hypotheses, indices0, indices1, 
#  var_e = 1, xmin = 0, xmax = 1, nugget = 1e-1, numCandidates = NULL, k = 4, p = 1, 
#  initN = NULL, initD = NULL, initD_indices = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
#  alpha_seq = 1, buffer_seq = 0, numDimsMax = NULL, 
#  genCandidates = 1, candidates = NULL, algorithm = 1, seed = NULL){
  
  
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
  # print(paste("finished ", 0, " out of ", numSeq, " steps", sep = ""))
  for(t in 1:numSeq){
    seed = seed + t
    set.seed(seed)
    Dt = add_MMEDgpvs(initD = D, y = y, type = type_hypotheses, l = l_hypotheses, 
                      indices0 = indices0, indices1 = indices1, var_e = var_e, 
                      N2 = N_seq[t], numCandidates = NULL, k = k, p = p, 
                      xmin = xmin, xmax = xmax, nugget = nugget, 
                      alpha = alpha_seq[t], buffer = buffer_seq[t], 
                      genCandidates = genCandidates, candidates = candidates, 
                      algorithm = algorithm)
    
    indt = Dt$indices
    yt = true_y[indt]
    
    # update D and y with new data
    D = rbind(D, Dt$addD)
    y = c(y, yt)
    indices = c(indices, indt)
    # print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
  }
  return(list("initD" = initD, 
              "addD" = D, 
              "D" = c(initD, D),
              "candidates" = candidates, 
              "indices" = indices))
}
