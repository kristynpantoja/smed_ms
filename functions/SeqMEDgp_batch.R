# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("covariance_functions.R")

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

obj_gp = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, y, error.var, type, l, 
  p = 1, k = 4, alpha = 1
){
  q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, error.var, type, l, p, 
            alpha)
  if(is.null(D)){ # when N2 = 1, and batch.idx != 1
    sum_q_D = sum(sapply(D, function(x_i) (1 / sqrt((x_i - candidate)^2))^k))
  } else{ # 
    sum_q_D = sum(sapply(D, function(x_i)
      (q_gp(x_i, Kinv0, Kinv1, initD, y, error.var, type, l, p,
            alpha) / sqrt((x_i - candidate)^2))^k))
  }
  result = q_cand^k * sum_q_D
  return(result)
}

# MMED_gp_batch, add_MED_ms_oneatatime_data_gp, add_MMEDgp_oneatatime
SeqMEDgp_batch = function(
  initD, y, type, l, error.var = 1, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, genCandidates = 1, 
  candidates = NULL, batch.idx = 1
){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # old_initD = initD
  
  # posterior distribution of beta
  if(is.null(nugget)){
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  } else{
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
  }
  
  initN = length(initD)
  ttlN = initN + N2
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, FUN = function(x) WNgp(
      x, Kinv0, Kinv1, initD, y, error.var, type, l))
    w_opt = which.max(w_candidates)
    xopt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))
  } else{
    is_x_max_in_initD = TRUE
  }
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) obj_gp(
        x, NULL, 
        Kinv0, Kinv1, initD, y, error.var, type, l, p, k, alpha))
    f_opt = which.min(f_min_candidates)
    
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = w_opt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = sapply(
        candidates, 
        function(x) obj_gp(
          x, D[1:(i - 1)], 
          Kinv0, Kinv1, initD, y, error.var, type, l, p, k, alpha))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt]
      # Update set of design points (D) and plot new point
      D[i] = xnew
      D_ind[i] = f_opt
    }
  }
  
  return(list(
    "initD" = initD, 
    "addD" = D, 
    "D" = c(initD, D),
    "candidates" = candidates, 
    "indices" = D_ind
  ))
}

########################################################################################################################################################################################

obj_gp.old = function(
  candidate, D, Kinv0, Kinv1, initD, y, error.var, type, l, p = 1, k = 4, alpha = 1
){
  result = q_gp(candidate, Kinv0, Kinv1, initD, y, error.var, type, l, p, 
                alpha)^k * 
  sum(sapply(D, function(x_i) (1 / sqrt((x_i - candidate)^2))^k))
  return(result)
}

# MMED_gp_batch, add_MED_ms_oneatatime_data_gp, add_MMEDgp_oneatatime
SeqMEDgp_batch.old = function(
  initD, y, type, l, error.var = 1, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, genCandidates = 1, 
  candidates = NULL, batch.idx = 1
){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # old_initD = initD
  
  # posterior distribution of beta
  if(is.null(nugget)){
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  } else{
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
  }
  
  # w_initD = sapply(initD, FUN = function(x) WNgp(
  #   x, Kinv0, Kinv1, initD, y, error.var, type, l))
  # # take out the values with wasserstein distance equal to 0
  # if(length(which(w_initD == 0)) != 0){
  #   initD = initD[-which(w_initD == 0)]
  #   y = y[-which(w_initD == 0)]
  #   w_initD = w_initD[-which(w_initD == 0)]
  #   initN = length(initD)
  #   # recalculate Kinv0 and Kinv1
  #   if(is.null(nugget)){
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  #   } else{
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
  #   }
  # }
  initN = length(initD)
  ttlN = initN + N2
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, FUN = function(x) WNgp(
      x, Kinv0, Kinv1, initD, y, error.var, type, l))
    w_opt = which.max(w_candidates)
    xopt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))
  } else{
    is_x_max_in_initD = TRUE
  }
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) obj_gp.old(
        x, initD, Kinv0, Kinv1, initD, y, error.var, type, l, p, k, alpha))
    f_opt = which.min(f_min_candidates)
    
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = w_opt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = sapply(candidates, function(x) obj_gp.old(x, D[1:(i - 1)], Kinv0, Kinv1, initD, y, error.var, type, l, p, k, alpha))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt]
      # Update set of design points (D) and plot new point
      D[i] = xnew
      D_ind[i] = f_opt
    }
  }
  
  return(list(
    "initD" = initD, 
    "addD" = D, 
    "D" = c(initD, D),
    "candidates" = candidates, 
    "indices" = D_ind
    ))
}






################################################
### SeqMED GP for Variable Selection, 1D vs 2D ###
################################################

# meant to be able to handle 2d dimensional input, for variable selection problem

obj_gpvs = function(
  candidate, D, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, error.var, 
  type, l, p = 1, k = 4, alpha = 1
){
  result = q_gpvs(
    candidate, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, error.var, 
    type, l, p, alpha)^k * 
    sum(apply(D, 1, function(x_i) (1 / sqrt(sum((x_i - candidate)^2)))^k))
  return(result)
}

#add_MMEDgpvs_oneatatime
SeqMEDgpvs_batch_oneatatime = function(
  initD, y, type = c(1, 1), l = c(0.1, 0.1), indices0, indices1, error.var = 1, 
  N2 = 11, numCandidates = 10^5, k = 4, p = 1, xmin = 0, xmax = 1, 
  nugget = NULL, alpha = NULL, genCandidates = 1, candidates = NULL, 
  batch.idx = 1
){
  if(typeof(initD) == "list") initD = as.matrix(initD)
  initN = dim(initD)[1]
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 
  # (in which case we don't want to use it, since 1/0 in q)
  old_initD = initD
  
  ttlN = initN + N2
  initD0 = NULL
  # posterior distribution of beta
  if(type[1] == type[2]){
    if(!is.null(indices0)){
      initD0 = initD[ , indices0, drop = FALSE]
    } else{
      warning("types are the same, but indices0 was not given; assumed to be V1")
      indices0 = 1
      initD0 = initD[ , indices0, drop = FALSE]
    }
    if(!is.null(indices1)){
      initD1 = initD[ , indices1, drop = FALSE]
    } else{
      warning("types are the same, but indices1 was not given; assumed to be all variables")
      indices1 = c(1:dim(initD)[2])
      initD0 = initD[ , indices1, drop = FALSE]
    }
    if(is.null(nugget)){
      Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]))
      Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]) + diag(rep(nugget, initN)))
      Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]) + diag(rep(nugget, initN)))
    }
  } else{
    if(is.null(nugget)){
      Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]))
      Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]) + diag(rep(nugget, initN)))
      Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]) + diag(rep(nugget, initN)))
    }
  }
  
  # -- Generate Candidate Points -- #
  if(!is.null(candidates)){
    numCandidates = length(candidates)
  } else{
    if(genCandidates == 1) {
      candidates_marginal = seq(from = xmin, to = xmax, length.out = numCandidates)
      candidates = expand.grid(candidates_marginal, candidates_marginal)
    }
    if(genCandidates == 2) {
      candidates_marginal = sort(runif(numCandidates, min = xmin, max = xmax))
      candidates = expand.grid(candidates_marginal, candidates_marginal)
    }
  }
  candidates = as.matrix(candidates)
  
  # -- Initialize 1st additional design point-- #
  D = matrix(rep(NA, N2 * 2), N2, 2)
  D_ind = rep(NA, N2)
  w_evals = apply(candidates, 1, FUN = function(x) WNgpvs(x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, error.var, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex, ]
  # I'm just gonna assume (and I think this is a safe assumption based on how Wasserstein was defined)
  # that the optimal point won't be one already in the initial design.
  D[1, ] = xopt
  D_ind[1] = xoptindex
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = apply(candidates, 1, function(x) obj_gpvs(x, D[1:(i - 1), , drop = FALSE], Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, error.var, type, l, p, k, alpha))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt, ]
      # Update set of design points (D) and plot new point
      D[i, ] = xnew
      D_ind[i] = f_opt
    }
  }
  
  return(list("initD" = initD, 
              "addD" = D, 
              "D" = c(initD, D),
              "candidates" = candidates, 
              "indices" = D_ind))
}

# add_MMEDgpvs
SeqMEDgpvs_batch = function(
  initD, y, type = c(1, 1), l = c(0.1, 0.1), indices0, indices1, error.var = 1, 
  N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, 
  genCandidates = 1, candidates = NULL, algorithm = 1, batch.idx = 1
){
  if(algorithm == 1){
    SeqMEDgpvs_batch_oneatatime(
      initD, y, type, l, indices0, indices1, error.var, N2, numCandidates, k, p, 
      xmin, xmax, nugget, alpha, genCandidates, candidates, batch.idx)
  } else{
    stop("invalid algorithm")
  }
}
