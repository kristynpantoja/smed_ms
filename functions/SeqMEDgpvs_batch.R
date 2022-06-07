#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

obj_newq_gpvs = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, initD0, initD1, y, 
  p = 1, k = 4, alpha = 1, buffer = 0, objective.type = 1, model0, model1
){
  ### objective.type == 1 adds a buffer to the wasserstein distance ###
  if(objective.type == 1 | objective.type %in% 
     c("buffer", "augdist", "augmenteddistance")){
    # q(x), x in C
    q_cand = q_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, 
                    y, p, alpha, buffer, model0, model1)
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(apply(initD, 1, function(x_i)
      (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, 
              alpha, buffer, model0, model1) / 
         sqrt(sum((x_i - candidate)^2)))^k))
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(apply(D, 1, function(x_i)  # no need for buffer here, either fyi
          (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, 
                  alpha, buffer, model0, model1) / 
             sqrt(sum((x_i - candidate)^2)))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 2 is just optimizing q(x) over x \in C ###
  if(objective.type == 2 | objective.type == "q"){
    # no need for buffer in this case
    q_cand = q_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, 
                    y, p, alpha, buffer = 0, model0, model1)
    if(is.null(D)){ # when N2 = 1, and batch.idx != 1
      result = q_cand^k
    } else{ # when N2 > 1
      sum_q_D = sum(apply(D, 1, function(x_i) # disregard initD
        (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, 
                alpha, buffer = 0, model0, model1) / 
           sqrt(sum((x_i - candidate)^2)))^k))
      result = q_cand^k * sum_q_D
    }
  }
  ### objective.type == 3 gives uniform charge to input points ###
  if(objective.type == 3 | objective.type %in% c("uniform", "spacefilling")){
    # q(x), x in C
    q_cand = q_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, 
                    y, p, alpha, buffer = 0, model0, model1)
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(apply(initD, 1, function(x_i) 
      (1 / sqrt(sum((x_i - candidate)^2)))^k)) # q = 1 for xi in this case
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(apply(D, 1, function(x_i) 
          (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, 
                  alpha, buffer = 0, model0, model1) / 
             sqrt(sum((x_i - candidate)^2)))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 4 caps q ###
  if(objective.type == 4 | objective.type %in% c("cap", "capq")){
    # q(x), x in C
    q_cand = qcap_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, y, p, alpha, 
                     model0, model1) # no capping needed
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    sum_q_D = sum(apply(initD, 1, function(x_i)
      (qcap_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha, 
               model0, model1) / 
         sqrt(sum((x_i - candidate)^2)))^k))
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(apply(D, 1, function(x_i)  # no capping needed here, either
          (qcap_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha, 
                   model0, model1) / 
             sqrt(sum((x_i - candidate)^2)))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### objective.type == 5 leaves out each x_i in D^{(t)} in q(x_i) calc ###
  if(objective.type == 5 | objective.type %in% 
     c("loo", "leaveout", "leaveoneout")){
    # q(x), x in C
    q_cand = q_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, y, p, alpha, 
                  buffer = 0,  model0, model1) # no need to leave anything out here, jsyk
    # the other terms in the summation
    # q(xi), xi in the observed set, D_t^c
    terms_q_D = rep(NA, nrow(initD))
    for(i in 1:nrow(initD)){ 
      terms_q_D[i] = (q_gpvs(
        initD[i, ], 
        Kinv0[-i, -i, drop = FALSE], Kinv1[-i, -i, drop = FALSE], 
        initD0[-i, , drop = FALSE], initD1[-i, , drop = FALSE], 
        y[-i], p, alpha, buffer = 0, model0, model1) / 
          sqrt(sum((initD[i, , drop = FALSE] - candidate)^2)))^k
    }
    sum_q_D = sum(terms_q_D)
    if(!is.null(D)){ # when N2 > 1
      # q(xi), xi in the unobserved design points D^{(t)}
      sum_q_D = sum_q_D + 
        sum(apply(D, 1, function(x_i) # no leave out necessary here, either
          (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha, buffer = 0, 
                model0, model1) / 
             sqrt(sum((x_i - candidate)^2)))^k))
    }
    result = q_cand^k * sum_q_D
  }
  ### no other objective.type options
  if(!(objective.type %in% c(1, 2, 3, 4, 5))){
    stop("obj_newq_gpvs: invalid objective.type value")
  }
  return(result)
}

SeqMEDgpvs_newq_batch = function(
  initD, y, N2 = 11, numCandidates = NULL, k = 4, p = 1, 
  xmin = 0, xmax = 1, alpha = NULL, candidates, 
  batch.idx = 1, buffer = 0, objective.type = 1, model0, model1
){
  
  if(!("matrix" %in% class(initD))) stop("SeqMEDgpvs_newq_batch: initD isn't a matrix!")
  initN = dim(initD)[1]
  if(length(y) != initN) stop("SeqMEDgpvs_newq_batch: length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary

  
  # select the variables in each model
  if(!is.null(model0$indices)){
    initD0 = initD[ , model0$indices, drop = FALSE]
  } else{
    stop("SeqMEDgpvs_keepq_batch: model0$indices is not given!")
  }
  if(!is.null(model1$indices)){
    initD1 = initD[ , model1$indices, drop = FALSE]
  } else{
    stop("SeqMEDgpvs_keepq_batch: model1$indices is not given!")
  }
  
  # posterior predictive distribution Kinv's
  if(is.null(model0$measurement.var)){
    Kinv0 = solve(getCov(
      X1 = initD0, X2 = initD0, type = model0$type, l = model0$l, p = model0$p, 
      signal.var = model0$signal.var))
  } else{
    Kinv0 = solve(getCov(
      X1 = initD0, X2 = initD0, type = model0$type, l = model0$l, p  = model0$p, 
      signal.var = model0$signal.var) + model0$measurement.var * diag(initN))
  }
  if(is.null(model1$measurement.var)){
    Kinv1 = solve(getCov(
      X1 = initD1, X2 = initD1, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var))
  } else{
    Kinv1 = solve(getCov(
      X1 = initD1, X2 = initD1, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var) + model1$measurement.var * diag(initN))
  }
  
  # get SeqMED
  D = matrix(rep(NA, N2 * 2), N2, 2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = apply(candidates, 1, FUN = function(x) WNgpvs(
      x, Kinv0, Kinv1, initD0, initD1, y, model0, model1))
    w_opt = which.max(w_candidates)
    x_w_opt = candidates[w_opt, , drop = FALSE]
    is_x_max_in_initD = any(apply(initD, 1, function(x) 
      all(as.numeric(x) == as.numeric(x_w_opt))))
  } else{
    is_x_max_in_initD = TRUE
  }
  
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = apply(
      candidates, 1, FUN = function(x) obj_newq_gpvs(
        x, NULL, Kinv0, Kinv1, initD, initD0, initD1, y, p, k, alpha, buffer, 
        objective.type, model0, model1))
    if(all(f_min_candidates == Inf)){
      stop("SeqMEDgpvs_newq_batch: all candidates result in objective function = Inf.")
    }
    f_opt = which.min(f_min_candidates)
    x_f_opt = candidates[f_opt, , drop = FALSE]
    # Update set of design points (D) and plot new point
    D[1, ] = x_f_opt
    D_ind[1] = f_opt
  } else{
    D[1, ] = x_w_opt
    D_ind[1] = w_opt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = apply(
        candidates, 1, 
        function(x) obj_newq_gpvs(
          x, D[1:(i - 1), , drop = FALSE], 
          Kinv0, Kinv1, initD, initD0, initD1, y, p, k, alpha, buffer, 
          objective.type, model0, model1))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt, , drop = FALSE]
      # Update set of design points (D) and plot new point
      D[i, ]  = xnew
      D_ind[i] = f_opt
    }
  }
  
  return(list(
    "initD" = initD, 
    "addD" = D, 
    "D" = rbind(initD, D),
    "candidates" = candidates, 
    "indices" = D_ind
  ))
}

obj_keepq_gpvs = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, initD0, initD1, y, 
  p = 1, k = 4, alpha = 1, buffer = 0, objective.type = 1, model0, model1, qs
){
  if(nrow(initD) != length(qs)) stop("obj_keepq_gpvs: nrow(initD) != length(qs)")
  # q(x), x in C
  
  q_cand = q_gpvs(candidate, Kinv0, Kinv1, initD0, initD1, 
                  y, p, alpha, buffer = 0, model0, model1)
  # the other terms in the summation
  # q(xi), xi in the observed set, D_t^c
  sum_q_D = sum((qs / sqrt(sum((initD - candidate)^2)))^k)
  # sum_q_D = sum(sapply(qs, function(q_i) 
  #   (q_i / sqrt((x_i - candidate)^2))^k)) # q = 1 for xi in this case #########
  if(!is.null(D)){ # when N2 > 1
    # q(xi), xi in the unobserved design points D^{(t)}
    sum_q_D = sum_q_D + 
      sum(apply(D, 1, function(x_i) 
        (q_gpvs(x_i, Kinv0, Kinv1, initD0, initD1, y, p, alpha, buffer = 0, 
              model0, model1) / 
           sqrt((x_i - candidate)^2))^k))
  }
  result = q_cand^k * sum_q_D
  
  return(data.frame(objectives = result, q.candidates = q_cand))
}

SeqMEDgpvs_keepq_batch = function(
  initD, y, N2 = 11, numCandidates = NULL, k = 4, p = 1, 
  xmin = 0, xmax = 1, alpha = NULL, candidates, 
  batch.idx = 1, buffer = 0, objective.type = 1, model0, model1, qs
){
  
  if(!("matrix" %in% class(initD))) stop("SeqMEDgpvs_keepq_batch: initD isn't a matrix!")
  initN = dim(initD)[1]
  if(length(y) != initN) stop("SeqMEDgpvs_keepq_batch: length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # select the variables in each model
  if(!is.null(model0$indices)){
    initD0 = initD[ , model0$indices, drop = FALSE]
  } else{
    stop("SeqMEDgpvs_keepq_batch: model0$indices is not given!")
  }
  if(!is.null(model1$indices)){
    initD1 = initD[ , model1$indices, drop = FALSE]
  } else{
    stop("SeqMEDgpvs_keepq_batch: model1$indices is not given!")
  }
  
  # posterior predictive distribution Kinv's
  if(is.null(model0$measurement.var)){
    Kinv0 = solve(getCov(
      X1 = initD0, X2 = initD0, type = model0$type, l = model0$l, p = model0$p, 
      signal.var = model0$signal.var))
  } else{
    Kinv0 = solve(getCov(
      X1 = initD0, X2 = initD0, type = model0$type, l = model0$l, p  = model0$p, 
      signal.var = model0$signal.var) + model0$measurement.var * diag(initN))
  }
  if(is.null(model1$measurement.var)){
    Kinv1 = solve(getCov(
      X1 = initD1, X2 = initD1, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var))
  } else{
    Kinv1 = solve(getCov(
      X1 = initD1, X2 = initD1, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var) + model1$measurement.var * diag(initN))
  }
  
  # get SeqMED
  D = matrix(rep(NA, N2 * 2), N2, 2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = apply(candidates, 1, FUN = function(x) WNgpvs(
      x, Kinv0, Kinv1, initD0, initD1, y, model0, model1))
    w_opt = which.max(w_candidates)
    x_w_opt = candidates[w_opt, , drop = FALSE]
    is_x_max_in_initD = any(apply(initD, 1, function(x) 
      all(as.numeric(x) == as.numeric(x_w_opt))))
  } else{
    is_x_max_in_initD = TRUE
  }
  
  # Find f_opt: minimum of f_min --- pulled this out of the below if statement
  f_min_candidates = apply(
    candidates, 1, FUN = function(x) obj_keepq_gpvs(
      x, NULL, Kinv0, Kinv1, initD, initD0, initD1, y, p, k, alpha, buffer, 
      objective.type, model0, model1, qs), simplify = FALSE) ####################################
  f_min_candidates = do.call("rbind", f_min_candidates)
  if(all(f_min_candidates$objectives == Inf)){
    stop("SeqMEDgp_batch: all candidates result in objective function = Inf.")
  }
  f_opt = which.min(f_min_candidates$objectives)
  x_f_opt = candidates[f_opt, , drop = FALSE]
  # Update set of design points (D) and plot new point
  if(is_x_max_in_initD){
    D[1, ] = x_f_opt
    D_ind[1] = f_opt
    q.new = f_min_candidates$q.candidates[f_opt]
  } else{
    D[1, ] = x_w_opt
    D_ind[1] = w_opt
    q.new = f_min_candidates$q.candidates[w_opt]
  }
  
  if(N2 > 1){
    if(N2 < 2){
      stop("SeqMEDgpvs_keepq_batch: To do a batch with N_batch > 2 using leave-out method, need to have at least 2 observations, to be able to leave one out!")
    }
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = apply(
        candidates, 1, FUN = function(x) obj_keepq_gpvs(
          x, D[1:(i - 1)], Kinv0, Kinv1, initD, initD0, initD1, y, p, k, alpha, 
          buffer, objective.type, model0, model1, qs), simplify = FALSE) ####################################
      f_min_candidates = do.call("rbind", f_min_candidates)
      f_opt = which.min(f_min_candidates$objectives)
      xnew = candidates[f_opt, , drop = FALSE]
      # Update set of design points (D) and plot new point
      D[i, ]  = xnew
      D_ind[i] = f_opt
      q.new = c(q.new, f_min_candidates$q.candidates[f_opt])
    }
  }
  
  return(list(
    "initD" = initD, 
    "addD" = D, 
    "D" = rbind(initD, D),
    "candidates" = candidates, 
    "indices" = D_ind,
    "q.new" = q.new
  ))
}
