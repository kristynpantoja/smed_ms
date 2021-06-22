# require("wasserstein_distance.R")
# require("charge_function_q.R")
# require("covariance_functions.R")

#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

obj_newq_gp = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, y, 
  p = 1, k = 4, alpha = 1, buffer = 0, objective.type = 1, model0, model1
){
  ### objective.type == 1 adds a buffer to the wasserstein distance ###
  if(objective.type == 1 | objective.type %in% 
     c("buffer", "augdist", "augmenteddistance")){
    # q(x), x in C
    q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, buffer, 
                  model0, model1) # no need for buffer here, jsyk
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(sapply(initD, function(x_i) 
      (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer, 
            model0, model1) / 
         sqrt((x_i - candidate)^2))^k))
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(sapply(D, function(x_i)  # no need for buffer here, either fyi
          (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer, 
                model0, model1) / 
             sqrt((x_i - candidate)^2))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 2 is just optimizing q(x) over x \in C ###
  if(objective.type == 2 | objective.type == "q"){
    # no need for buffer in this case
    q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                  model0, model1)
    if(is.null(D)){ # when N2 = 1, and batch.idx != 1
      result = q_cand^k
    } else{ # when N2 > 1
      sum_q_D = sum(sapply(D, function(x_i) # disregard initD
        (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
              model0, model1) / 
           sqrt((x_i - candidate)^2))^k))
      result = q_cand^k * sum_q_D
    }
  }
  ### objective.type == 3 gives uniform charge to input points ###
  if(objective.type == 3 | objective.type %in% c("uniform", "spacefilling")){
    # q(x), x in C
    q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                  model0, model1)
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(sapply(initD, function(x_i) 
      (1 / sqrt((x_i - candidate)^2))^k)) # q = 1 for xi in this case
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(sapply(D, function(x_i) 
          (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                model0, model1) / 
             sqrt((x_i - candidate)^2))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 4 caps q ###
  if(objective.type == 4 | objective.type %in% c("cap", "capq")){
    # q(x), x in C
    q_cand = qcap_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, 
                     model0, model1) # no capping needed
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(sapply(initD, function(x_i) 
      (qcap_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, 
               model0, model1) / 
         sqrt((x_i - candidate)^2))^k))
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(sapply(D, function(x_i)  # no capping needed here, either
          (qcap_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, 
                   model0, model1) / 
             sqrt((x_i - candidate)^2))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 5 leaves out each x_i in D^{(t)} in q(x_i) calc ###
  if(objective.type == 5 | objective.type %in% 
     c("loo", "leaveout", "leaveoneout")){
    # q(x), x in C
    q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                  model0, model1) # no need to leave anything out here, jsyk
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    terms_q_D = rep(NA, length(initD))
    for(i in 1:length(initD)){
      terms_q_D[i] = (q_gp(
        initD[i], Kinv0[-i, -i, drop = FALSE], Kinv1[-i, -i, drop = FALSE], 
        initD[-i], y[-i], p, alpha, buffer = 0, model0, model1) / 
          sqrt((initD[i] - candidate)^2))^k
    }
    sum_q_D = sum(terms_q_D)
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(sapply(D, function(x_i) # no leave out necessary here, either
          (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                model0, model1) / 
             sqrt((x_i - candidate)^2))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### no other objective.type options
  if(!(objective.type %in% c(1, 2, 3, 4, 5))){
    stop("obj_newq_gp: invalid objective.type value")
  }
  return(result)
}

SeqMEDgp_newq_batch = function(
  initD, y, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, alpha = NULL, candidates = NULL, 
  batch.idx = 1, buffer = 0, objective.type = 1, model0, model1
){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # old_initD = initD
  
  # posterior distribution of beta
  if(is.null(model0$measurement.var)){
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l))
  } else{
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l) + 
                    sqrt(model0$measurement.var) * diag(initN))
  }
  if(is.null(model1$measurement.var)){
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l))
  } else{
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l) + 
                    sqrt(model1$measurement.var) * diag(initN))
  }
  
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, FUN = function(x) WNgp(
      x, Kinv0, Kinv1, initD, y, model0, model1))
    w_opt = which.max(w_candidates)
    x_w_opt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == x_w_opt))
  } else{
    is_x_max_in_initD = TRUE
  }
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) obj_newq_gp(
        x, NULL, 
        Kinv0, Kinv1, initD, y, p, k, alpha, buffer, objective.type, 
        model0, model1))
    if(all(f_min_candidates == Inf)){
      stop("SeqMEDgp_batch: all candidates result in objective function = Inf.")
    }
    f_opt = which.min(f_min_candidates)
    x_f_opt = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = x_f_opt
    D_ind[1] = f_opt
  } else{
    D[1] = x_w_opt
    D_ind[1] = w_opt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = sapply(
        candidates, 
        function(x) obj_newq_gp(
          x, D[1:(i - 1)], 
          Kinv0, Kinv1, initD, y, p, k, alpha, buffer, objective.type, 
          model0, model1))
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

obj_keepq_gp = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, y, 
  p = 1, k = 4, alpha = 1, buffer = 0, objective.type = 1, model0, model1, qs
){
  if(length(initD) != length(qs)) stop("obj_keepq_gp: length(initD) != length(qs)")
    # q(x), x in C
    q_cand = q_gp(candidate, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                  model0, model1)
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum((qs / sqrt((initD - candidate)^2))^k)
    # sum_q_D = sum(sapply(qs, function(q_i) 
    #   (q_i / sqrt((x_i - candidate)^2))^k)) # q = 1 for xi in this case #########
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(sapply(D, function(x_i) 
          (q_gp(x_i, Kinv0, Kinv1, initD, y, p, alpha, buffer = 0, 
                model0, model1) / 
             sqrt((x_i - candidate)^2))^k))
    }
    result = q_cand^k * sum_q_D
  
  return(data.frame(objectives = result, q.candidates = q_cand))
}

SeqMEDgp_keepq_batch = function(
  initD, y, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
  xmin = 0, xmax = 1, alpha = NULL, candidates = NULL, 
  batch.idx = 1, buffer = 0, objective.type = 1, model0, model1, qs
){
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  if(is.null(qs)) stop("SeqMEDgp_keepq_batch: qs == NULL!")
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # old_initD = initD
  
  # posterior distribution of beta
  if(is.null(model0$measurement.var)){
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l))
  } else{
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l) + 
                    sqrt(model0$measurement.var) * diag(initN))
  }
  if(is.null(model1$measurement.var)){
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l))
  } else{
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l) + 
                    sqrt(model1$measurement.var) * diag(initN))
  }
  
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, FUN = function(x) WNgp(
      x, Kinv0, Kinv1, initD, y, model0, model1))
    w_opt = which.max(w_candidates)
    x_w_opt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == x_w_opt))
  } else{
    is_x_max_in_initD = TRUE
  }
  # Find f_opt: minimum of f_min --- pulled this out of the below if statement
  f_min_candidates = lapply(
    candidates, 
    function(x) obj_keepq_gp(
      x, NULL, 
      Kinv0, Kinv1, initD, y, p, k, alpha, buffer, objective.type, 
      model0, model1, qs))
  f_min_candidates = do.call("rbind", f_min_candidates)
  if(all(f_min_candidates$objectives == Inf)){
    stop("SeqMEDgp_batch: all candidates result in objective function = Inf.")
  }
  f_opt = which.min(f_min_candidates$objectives)
  x_f_opt = candidates[f_opt]
  # Update set of design points (D) and plot new point
  if(is_x_max_in_initD){
    D[1] = x_f_opt
    D_ind[1] = f_opt
    q.new = f_min_candidates$q.candidates[f_opt]
  } else{
    D[1] = x_w_opt
    D_ind[1] = w_opt
    q.new = f_min_candidates$q.candidates[w_opt]
  }
  
  if(N2 > 1){
    if(N2 < 2){
      stop("SeqMED_keepq_batch: To do a batch with N_batch > 2 using leave-out method, need to have at least 2 observations, to be able to leave one out!")
    }
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = lapply(
        candidates, 
        function(x) obj_keepq_gp(
          x, D[1:(i - 1)], 
          Kinv0, Kinv1, initD, y, p, k, alpha, buffer, objective.type, 
          model0, model1, qs))
      f_min_candidates = do.call("rbind", f_min_candidates)
      f_opt = which.min(f_min_candidates$objectives)
      xnew = candidates[f_opt]
      # Update set of design points, D
      D[i] = xnew
      D_ind[i] = f_opt
      q.new = c(q.new, f_min_candidates$q.candidates[f_opt])
    }
  }
  
  return(list(
    "initD" = initD, 
    "addD" = D, 
    "D" = c(initD, D),
    "candidates" = candidates, 
    "indices" = D_ind,
    "q.new" = q.new
  ))
}

########################################################################################################################################################################################

# obj_gp.old = function(
#   candidate, D, Kinv0, Kinv1, initD, y, signal.var, type, l, p = 1, k = 4, alpha = 1
# ){
#   result = q_gp(candidate, Kinv0, Kinv1, initD, y, signal.var, type, l, p, 
#                 alpha)^k * 
#     sum(sapply(D, function(x_i) (1 / sqrt((x_i - candidate)^2))^k))
#   return(result)
# }
# 
# # MMED_gp_batch, add_MED_ms_oneatatime_data_gp, add_MMEDgp_oneatatime
# SeqMEDgp_batch.old = function(
#   initD, y, type, l, signal.var = 1, N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
#   xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, genCandidates = 1, 
#   candidates = NULL, batch.idx = 1
# ){
#   initN = length(initD)
#   if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
#   
#   # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
#   # turns out, that's not necessary
#   
#   # old_initD = initD
#   
#   # posterior distribution of beta
#   if(is.null(nugget)){
#     Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
#     Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
#   } else{
#     Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
#     Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
#   }
#   
#   # w_initD = sapply(initD, FUN = function(x) WNgp(
#   #   x, Kinv0, Kinv1, initD, y, signal.var, type, l))
#   # # take out the values with wasserstein distance equal to 0
#   # if(length(which(w_initD == 0)) != 0){
#   #   initD = initD[-which(w_initD == 0)]
#   #   y = y[-which(w_initD == 0)]
#   #   w_initD = w_initD[-which(w_initD == 0)]
#   #   initN = length(initD)
#   #   # recalculate Kinv0 and Kinv1
#   #   if(is.null(nugget)){
#   #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]))
#   #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]))
#   #   } else{
#   #     Kinv0 = solve(getCov(initD, initD, type[1], l[1]) + diag(rep(nugget, initN)))
#   #     Kinv1 = solve(getCov(initD, initD, type[2], l[2]) + diag(rep(nugget, initN)))
#   #   }
#   # }
#   initN = length(initD)
#   ttlN = initN + N2
#   
#   # -- Generate Candidate Points -- #
#   if(is.null(candidates)){
#     if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
#     if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
#   }
#   
#   D = rep(NA, N2)
#   D_ind = rep(NA, N2)
#   if(batch.idx == 1){
#     # -- Initialize 1st additional design point-- #
#     w_candidates = sapply(candidates, FUN = function(x) WNgp(
#       x, Kinv0, Kinv1, initD, y, signal.var, type, l))
#     w_opt = which.max(w_candidates)
#     xopt = candidates[w_opt]
#     is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))
#   } else{
#     is_x_max_in_initD = TRUE
#   }
#   if(is_x_max_in_initD){
#     # Find f_opt: minimum of f_min
#     f_min_candidates = sapply(
#       candidates, 
#       function(x) obj_gp.old(
#         x, initD, Kinv0, Kinv1, initD, y, signal.var, type, l, p, k, alpha))
#     f_opt = which.min(f_min_candidates)
#     
#     xnew = candidates[f_opt]
#     # Update set of design points (D) and plot new point
#     D[1] = xnew
#     D_ind[1] = f_opt
#   } else{
#     D[1] = xopt
#     D_ind[1] = w_opt
#   }
#   
#   if(N2 > 1){
#     for(i in 2:N2){
#       # Find f_opt: minimum of f_min
#       f_min_candidates = sapply(candidates, function(x) obj_gp.old(x, D[1:(i - 1)], Kinv0, Kinv1, initD, y, signal.var, type, l, p, k, alpha))
#       f_opt = which.min(f_min_candidates)
#       xnew = candidates[f_opt]
#       # Update set of design points (D) and plot new point
#       D[i] = xnew
#       D_ind[i] = f_opt
#     }
#   }
#   
#   return(list(
#     "initD" = initD, 
#     "addD" = D, 
#     "D" = c(initD, D),
#     "candidates" = candidates, 
#     "indices" = D_ind
#   ))
# }






# ##################################################
# ### SeqMED GP for Variable Selection, 1D vs 2D ###
# ##################################################
# 
# # meant to be able to handle 2d dimensional input, for variable selection problem
# 
# obj_gpvs = function(
#   candidate, D, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, 
#   type, l, p = 1, k = 4, alpha = 1
# ){
#   result = q_gpvs(
#     candidate, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, 
#     type, l, p, alpha)^k * 
#     sum(apply(D, 1, function(x_i) (1 / sqrt(sum((x_i - candidate)^2)))^k))
#   return(result)
# }
# 
# #add_MMEDgpvs_oneatatime
# SeqMEDgpvs_batch_oneatatime = function(
#   initD, y, type = c(1, 1), l = c(0.1, 0.1), indices0, indices1, signal.var = 1, 
#   N2 = 11, numCandidates = 10^5, k = 4, p = 1, xmin = 0, xmax = 1, 
#   nugget = NULL, alpha = NULL, genCandidates = 1, candidates = NULL, 
#   batch.idx = 1
# ){
#   if(typeof(initD) == "list") initD = as.matrix(initD)
#   initN = dim(initD)[1]
#   if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
#   
#   # check if any points in initD give Wasserstein distance of 0 
#   # (in which case we don't want to use it, since 1/0 in q)
#   old_initD = initD
#   
#   ttlN = initN + N2
#   initD0 = NULL
#   # posterior distribution of beta
#   if(type[1] == type[2]){
#     if(!is.null(indices0)){
#       initD0 = initD[ , indices0, drop = FALSE]
#     } else{
#       warning("types are the same, but indices0 was not given; assumed to be V1")
#       indices0 = 1
#       initD0 = initD[ , indices0, drop = FALSE]
#     }
#     if(!is.null(indices1)){
#       initD1 = initD[ , indices1, drop = FALSE]
#     } else{
#       warning("types are the same, but indices1 was not given; assumed to be all variables")
#       indices1 = c(1:dim(initD)[2])
#       initD0 = initD[ , indices1, drop = FALSE]
#     }
#     if(is.null(nugget)){
#       Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]))
#       Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]))
#     } else{
#       Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]) + diag(rep(nugget, initN)))
#       Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]) + diag(rep(nugget, initN)))
#     }
#   } else{
#     if(is.null(nugget)){
#       Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]))
#       Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]))
#     } else{
#       Kinv0 = solve(getCov(initD0, initD0, type[1], l[1]) + diag(rep(nugget, initN)))
#       Kinv1 = solve(getCov(initD1, initD1, type[2], l[2]) + diag(rep(nugget, initN)))
#     }
#   }
#   
#   # -- Generate Candidate Points -- #
#   if(!is.null(candidates)){
#     numCandidates = length(candidates)
#   } else{
#     if(genCandidates == 1) {
#       candidates_marginal = seq(from = xmin, to = xmax, length.out = numCandidates)
#       candidates = expand.grid(candidates_marginal, candidates_marginal)
#     }
#     if(genCandidates == 2) {
#       candidates_marginal = sort(runif(numCandidates, min = xmin, max = xmax))
#       candidates = expand.grid(candidates_marginal, candidates_marginal)
#     }
#   }
#   candidates = as.matrix(candidates)
#   
#   # -- Initialize 1st additional design point-- #
#   D = matrix(rep(NA, N2 * 2), N2, 2)
#   D_ind = rep(NA, N2)
#   w_evals = apply(candidates, 1, FUN = function(x) WNgpvs(x, Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, type, l))
#   xoptindex = which.max(w_evals)
#   xopt = candidates[xoptindex, ]
#   # I'm just gonna assume (and I think this is a safe assumption based on how Wasserstein was defined)
#   # that the optimal point won't be one already in the initial design.
#   D[1, ] = xopt
#   D_ind[1] = xoptindex
#   
#   if(N2 > 1){
#     for(i in 2:N2){
#       # Find f_opt: minimum of f_min
#       f_min_candidates = apply(candidates, 1, function(x) obj_gpvs(x, D[1:(i - 1), , drop = FALSE], Kinv0, Kinv1, indices0, indices1, initD0, initD1, y, signal.var, type, l, p, k, alpha))
#       f_opt = which.min(f_min_candidates)
#       xnew = candidates[f_opt, ]
#       # Update set of design points (D) and plot new point
#       D[i, ] = xnew
#       D_ind[i] = f_opt
#     }
#   }
#   
#   return(list("initD" = initD, 
#               "addD" = D, 
#               "D" = c(initD, D),
#               "candidates" = candidates, 
#               "indices" = D_ind))
# }
# 
# # add_MMEDgpvs
# SeqMEDgpvs_batch = function(
#   initD, y, type = c(1, 1), l = c(0.1, 0.1), indices0, indices1, signal.var = 1, 
#   N2 = 11, numCandidates = 10^5, k = 4, p = 1, 
#   xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, 
#   genCandidates = 1, candidates = NULL, algorithm = 1, batch.idx = 1
# ){
#   if(algorithm == 1){
#     SeqMEDgpvs_batch_oneatatime(
#       initD, y, type, l, indices0, indices1, signal.var, N2, numCandidates, k, p, 
#       xmin, xmax, nugget, alpha, genCandidates, candidates, batch.idx)
#   } else{
#     stop("invalid algorithm")
#   }
# }
