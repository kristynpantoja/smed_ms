#################################################
### SeqMED GP, one-at-a-time greedy algorithm ###
#################################################

obj_newq_gp = function(
  candidate, D = NULL, Kinv0, Kinv1, initD, y, 
  p = 1, k = NULL, alpha = 1, model0, model1
){
  if(is.null(k)) k = 4 * p
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
  return(result)
}

SeqMEDgp_newq_batch = function(
  initD, y, N2 = 1, numCandidates = NULL, k = NULL, p = 1, 
  xmin = -1, xmax = 1, alpha = NULL, candidates = NULL, 
  batch.idx = 1, model0, model1
){
  if(is.null(k)) k = 4 * p
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # posterior predictive distribution Kinv's
  if(is.null(model0$measurement.var)){
    Kinv0 = solve(getCov(
      X1 = initD, X2 = initD, type = model0$type, l = model0$l, p = model0$p, 
      signal.var = model0$signal.var))
  } else{
    Kinv0 = solve(getCov(
      X1 = initD, X2 = initD, type = model0$type, l = model0$l, p  = model0$p, 
      signal.var = model0$signal.var) + model0$measurement.var * diag(initN))
  }
  if(is.null(model1$measurement.var)){
    Kinv1 = solve(getCov(
      X1 = initD, X2 = initD, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var))
  } else{
    Kinv1 = solve(getCov(
      X1 = initD, X2 = initD, type = model1$type, l = model1$l, p = model1$p, 
      signal.var = model1$signal.var) + model1$measurement.var * diag(initN))
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
        x, NULL, Kinv0, Kinv1, initD, y, p, k, alpha, model0, model1))
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
          x, D[1:(i - 1)], Kinv0, Kinv1, initD, y, p, k, alpha, model0, model1))
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
