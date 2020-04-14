# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("covariance_functions.R")

###############################################
### MMED GP, one-at-a-time greedy algorithm ###
###############################################

# old name: f_min_data_gp
f_min_gp = function(candidate, D, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer){
  result = q_gp(candidate, Kinv0, Kinv1, initD, y, var_e, type, l, p, 
                     alpha, buffer)^k * 
    sum(sapply(D, function(x_i) (q_gp(x_i, Kinv0, Kinv1, initD, y, var_e, type, l, p, 
                                           alpha, buffer) / sqrt((x_i - candidate)^2))^k))
  return(result)
}

# old name: add_MED_ms_oneatatime_data_gp
add_MMEDgp_oneatatime = function(initD, y, type, l, var_e, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
                                         xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                         genCandidates = 1, candidates = NULL){
  if(is.null(dim(initD))){
    initN = length(initD)
  } else{
    initN = dim(initD)[1]
  }
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  
  # posterior distribution of beta
  if(is.null(nugget)){
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  } else{
    Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
    Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
  }
  
  # check for any bad initial points in initD
  # w_initD = sapply(initD, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  # if(length(which(w_initD == 0)) != 0){
  #   initD = initD[-which(w_initD == 0)]
  #   y = y[-which(w_initD == 0)]
  #   w_initD = w_initD[-which(w_initD == 0)]
  #   if(is.null(nugget)){
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
  #   } else{
  #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, length(initD))))
  #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, length(initD))))
  #   }
  # }
  
  ttlN = initN + N2
  
  
  # -- Generate Candidate Points -- #
  if(!is.null(candidates)){
    numCandidates = length(candidates)
  } else{
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  # remove candidate points that are already in the design
  # this is necessary bc they lead to a posterior predictive distribution with
  # standard deviation/covariance = 0 and means equal to a design point
  # -> W = 0 -> division by 0 in q evaluation
  # candidates = candidates[!(candidates %in% x_train)] # actually not necessray, 
  # considering our edit Wasserstein_distance_gp = 0 when can't take sqrt
  
  # -- Initialize 1st additional design point-- #
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  w_evals = sapply(candidates, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, initD, y, var_e, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex]
  is_xopt_in_initD = any(sapply(initD, function(x) x == xopt)) # give tolerance?
  if(is_xopt_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_gp(x, initD, Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = xoptindex
  }
  
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(candidates, function(x) f_min_gp(x, D[1:(i - 1)], Kinv0, Kinv1, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i] = xnew
    D_ind[i] = f_opt
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD, 
              "candidates" = candidates, "indices" = D_ind))
}






################################################
### MMED GP for Variable Selection, 1D vs 2D ###
################################################

# meant to be able to handle 2d dimensional input, for variable selection problem
f_min_gpvs = function(candidate, D, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l, p, alpha, buffer){
  result = q_gpvs(candidate, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l, p, 
                      alpha, buffer)^k * 
    sum(apply(D, 1, function(x_i) (q_gpvs(x_i, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l, p, 
                                            alpha, buffer) / sqrt(sum((x_i - candidate)^2)))^k))
  return(result)
}

# old name: add_MED_ms_oneatatime_data_gp2
add_MMEDgpvs_oneatatime = function(initD, y, type = c(1, 1), l = c(0.1, 0.1), subdim = NULL, var_e = 1, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
                                   xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                   genCandidates = 1, candidates = NULL){
  # later, put some checks that things (like expand.grid) are matrices, NOT data frames!
  initD = unname(as.matrix(initD))
  initN = dim(initD)[1]
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  old_initD = initD
  ttlN = initN + N2
  subinitD = NULL
  # posterior distribution of beta
  if(type[1] == type[2]){
    if(!is.null(subdim)){
      subinitD = initD[ , subdim, drop = FALSE]
    } else{
      warning("types are the same, but subdim was not given; assumed to be V1")
      subinitD = initD[ , 1, drop = FALSE]
    }
    if(is.null(nugget)){
      Kinv0 = solve(getCov(subinitD, subinitD, type[1], l[1]))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(subinitD, subinitD, type[1], l[1]) + diag(rep(nugget, dim(subinitD)[1])))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, dim(initD)[1])))
    }
  } else{
    if(is.null(nugget)){
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, dim(initD)[1])))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, dim(initD)[1])))
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
  candidates = unname(as.matrix(candidates, dimnames = NULL))
  
  # -- Initialize 1st additional design point-- #
  D = matrix(rep(NA, N2 * 2), N2, 2)
  D_ind = rep(NA, N2)
  w_evals = apply(candidates, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex, ]
  # I'm just gonna assume (and I think this is a safe assumption based on how Wasserstein was defined)
  # that the optimal point won't be one already in the initial design.
  D[1, ] = xopt
  D_ind[1] = xoptindex
  
  if(N2 > 1){
    # temporarily, since uses global variables #########################################
    # quilt.plot(x_grid, y_seq)
    # plot(x_grid[ , 1], y_seq, ylim = range(y_seq))
    # order_truef = order(x_grid[,1])
    # lines(x_grid[order_truef,1], y_seq_mat_grid[order_truef,sim_index])
    # points(D[1, 1], y_seq[D_ind[1]], col = 2)
    #####################################################################################
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = apply(candidates, 1, function(x) f_min_gpvs(x, D[1:(i - 1), , drop = FALSE], Kinv0, Kinv1, subdim, subinitD, initD, y, var_e, type, l, p, alpha, buffer))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt, ]
      # Update set of design points (D) and plot new point
      D[i, ] = xnew
      D_ind[i] = f_opt
      # print(D[i, ])
      # temporarily, since uses global variables ##########################################
      # points(D[i, 1], y_seq[D_ind[i]], col = 2)
      #####################################################################################
    }
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD, 
              "candidates" = candidates, "indices" = D_ind))
}




################################################
### MMED GP for Variable Selection, all dims ###
################################################

add_MED_ms_oneatatime_data_gp3 = function(initD, y, type, l, subdim = NULL, var_e, N2 = 11, numCandidates = 10^5, k = 4, p = 2, 
                                          xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                          genCandidates = 1, candidates = NULL){
  # subdim is the vector of dimensions in the subspace
  
  # later, put some checks that things (like expand.grid) are matrices, NOT data frames!
  
  initN = dim(initD)[1]
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  
  old_initD = initD
  ttlN = initN + N2
  subinitD = NULL
  # posterior distribution of beta
  if(type[1] == type[2]){
    if(!is.null(subdim)){
      subinitD = initD[ , subdim, drop = FALSE]
    } else{
      warning("types are the same, but subdim was not given; assumed to be V1")
      subinitD = initD[ , 1, drop = FALSE]
    }
    if(is.null(nugget)){
      Kinv0 = solve(getCov(subinitD, subinitD, type[1], l[1]))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(subinitD, subinitD, type[1], l[1]) + diag(rep(nugget, dim(subinitD)[1])))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, dim(initD)[1])))
    }
  } else{
    if(is.null(nugget)){
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
    } else{
      Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, dim(initD)[1])))
      Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, dim(initD)[1])))
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
  w_evals = apply(candidates, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subinitD, initD, y, var_e, type, l))
  xoptindex = which.max(w_evals)
  xopt = candidates[xoptindex, ]
  # I'm just gonna assume (and I think this is a safe assumption based on how Wasserstein was defined)
  # that the optimal point won't be one already in the initial design.
  D[1, ] = xopt
  D_ind[1] = xoptindex
  
  
  # temporarily, since uses global variables #########################################
  # quilt.plot(x_grid, y_seq)
  # plot(x_grid[ , 1], y_seq, ylim = range(y_seq))
  # order_truef = order(x_grid[,1])
  # lines(x_grid[order_truef,1], y_seq_mat_grid[order_truef,sim_index])
  # points(D[1, 1], y_seq[D_ind[1]], col = 2)
  #####################################################################################
  for(i in 2:N2){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_gpvs(x, D[1:(i - 1), , drop = FALSE], Kinv0, Kinv1, subinitD, initD, y, var_e, type, l, p, alpha, buffer))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[i, ] = xnew
    D_ind[i] = f_opt
    # print(D[i, ])
    # temporarily, since uses global variables ##########################################
    # points(D[i, 1], y_seq[D_ind[i]], col = 2)
    #####################################################################################
  }
  
  return(list("initD" = old_initD, "addD" = D, "updatedD" = c(old_initD, D), "q_initD" = initD, 
              "candidates" = candidates, "indices" = D_ind))
}
