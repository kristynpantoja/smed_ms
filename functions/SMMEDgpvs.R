generate_SMMEDgpvs_oneatatime = function(type_true = NULL, l_true = NULL, subdim_true = NULL, true_y = NULL, 
                                         type_hypotheses, l_hypotheses, subdim_hypothesis, 
                                         var_e = 1, xmin = 0, xmax = 1, nugget = 1e-1, numCandidates = NULL, k = 4, p = 1, 
                                         initN = NULL, initD = NULL, initD_indices = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
                                         alpha_seq = 1, buffer_seq = 0, candidates = NULL, numDimsMax = NULL, genCandidates = 1, 
                                         seed = NULL){
  initD = as.matrix(initD)
  if(is.null(candidates) & is.null(true_y)){
    if(is.null(subdim_true)) subdim_true = 1:numDimsMax
    if(length(subdim_true) == 1) if(subdim_true != 1) subdim_true = 1:subdim_true
    if(length(subdim_hypothesis) == 1) if(subdim_hypothesis != 1) subdim_hypothesis = 1:subdim_hypothesis
    if(max(subdim_true) <= max(subdim_hypothesis)) stop("subdim_true is either less than or equal to the subsetting hypothesis")
    # get candidates, or discretized space
    if(is.null(candidates)){
      if(is.null(numDimsMax)){
        xseq = seq(from = xmin, to = xmax, length.out = numCandidates)
        xseqlist = list()
        for(i in 1:numDimsMax) xseqlist[[i]] = xseq 
        candidates = as.matrix(expand.grid(xseqlist))
      }
    }
  }
  candidates = as.matrix(candidates)
  # sequence checks
  if(length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
  if(length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(!is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  # check initD, or make it
  if(is.null(initD)){
    if(is.null(initN)){
      warning("initN and initD not given; going to generate N_seq[1] initial points")
      initD_indices = sample(1:dim(candidates)[1], N_seq[1])
      initD = candidates[initD_indices, ]
    } else{
      initD_indices = sample(1:dim(candidates)[1], initN)
      initD = candidates[initD_indices, ]
    }
  }
  # get true_y if not given (but hopefully it is)
  if(is.null(true_y)){
    null_cov = getCov(candidates[ , subdim_true], candidates[ , subdim_true], type_true, l_true)
    null_mean = rep(0, numCandidates^length(subdim_true))
    true_y = as.vector(rmvnorm(n = 1, mean = null_mean, sigma = null_cov))
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
  indices = initD_indices
  print(paste("finished ", 0, " out of ", numSeq, " steps", sep = ""))
  for(t in 1:numSeq){
    seed = seed + t
    set.seed(seed)
    Dt = add_MMEDgpvs_oneatatime(initD = D, y = y, type = type_hypotheses, l = l_hypotheses, subdim = subdim_hypothesis, 
                                 var_e = var_e, N2 = N_seq[t], numCandidates = NULL, k = k, p = p, 
                                 xmin = xmin, xmax = xmax, nugget = nugget, alpha = alpha_seq[t], buffer = buffer_seq[t], 
                                 genCandidates = genCandidates, candidates = candidates)
    ###
    initD = D
    y
    type = type_hypotheses
    l = l_hypotheses
    subdim = subdim_hypothesis
    var_e = 1
    N2 = N_seq[t]
    numCandidates = NULL
    k
    p
    xmin
    xmax
    nugget
    alpha = alpha_seq[t]
    buffer = buffer_seq[t]
    genCandidates
    candidates
    ###
    
    indt = Dt$indices
    yt = true_y[indt]
    
    # update D and y with new data
    D = rbind(D, Dt$addD)
    y = c(y, yt)
    indices = rbind(indices, indt)
    print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
  }
  return(list("D" = D, 
              "indices" = indices,
              "y" = y,
              "candidates" = candidates))
}

generate_SMMEDgpvs = function(type_true = NULL, l_true = NULL, subdim_true = NULL, true_y = NULL, 
                              type_hypotheses, l_hypotheses, subdim_hypothesis, 
                              var_e = 1, xmin = 0, xmax = 1, nugget = 1e-1, numCandidates = NULL, k = 4, p = 1, 
                              initN = NULL, initD = NULL, initD_indices = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
                              alpha_seq = 1, buffer_seq = 0, candidates = NULL, numDimsMax = NULL, genCandidates = 1, 
                              seed = NULL, algorithm = 1){
  if(algorithm == 1){
    generate_SMMEDgpvs_oneatatime(type_true, l_true, subdim_true, true_y, 
                                  type_hypotheses, l_hypotheses, subdim_hypothesis, 
                                  var_e, xmin, xmax, nugget, numCandidates, k, p, 
                                  initN, initD, initD_indices, inity, numSeq, N_seq, 
                                  alpha_seq, buffer_seq, candidates, numDimsMax, genCandidates, 
                                  seed)
  }
}
