# objective function 
f_min_vs = function( 
  candidate, D, model0, model1, postmean0, postmean1, postvar0, postvar1, 
  error.var, p, k, alpha
){
  result = q_vs(
    candidate, model0, model1, postmean0, postmean1, postvar0, postvar1, 
    error.var, p, alpha)^k * 
    sum(apply(D, 1, function(x_i) 
      (q_vs(x_i, model0, model1, postmean0, postmean1, postvar0, postvar1, 
            error.var, p, alpha) / sqrt(sum((x_i - candidate)^2)))^k))
  return(result)
}

# batch of seqN new points
SeqMEDvs_batch = function(
  initD, inity, model0, model1, error.var, N2 = 11, candidates, 
  dimX, xmin = -1, xmax = 1, k = 4, p = 1, alpha = 1, batch.idx
){
  initN = nrow(initD)
  if(length(inity) != initN) stop("length of inity does not match length of initial input data, initD")
  
  # posterior distributions of beta
  postbeta0 = getBetaPosterior(
    y = inity, X = initD[ , indices0, drop = FALSE], 
    beta.mean = model0$beta.mean, beta.var = model0$beta.var, 
    error.var = error.var)
  postbeta1 = getBetaPosterior(
    y = inity, X = initD[ , indices1, drop = FALSE], 
    beta.mean = model1$beta.mean, beta.var = model1$beta.var, 
    error.var = error.var)
  postvar0 = postbeta0$var
  postmean0 = postbeta0$mean
  postvar1 = postbeta1$var
  postmean1 = postbeta1$mean
  
  # check if any points in initD give Wasserstein distance of 0
  #   if there are any such points, remove them so that TPE does not explode
  #   (since 1/0 in q)
  w_initD = apply(initD, 1, FUN = function(x) WNlmvs(
      x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var))
  if(length(which(w_initD == 0)) != 0){
    initD = initD[-which(w_initD == 0), , drop = FALSE]
    y = inity[-which(w_initD == 0)]
    w_initD = w_initD[-which(w_initD == 0)]
    # recalculate posterior distributions of beta
    postbeta0 = getBetaPosterior(
      y = inity, X = initD[ , indices0], model0$beta.mean, model0$beta.var, 
      error.var)
    postbeta1 = getBetaPosterior(
      y = inity, X = initD[ , indices1], model1$beta.mean, model1$beta.var, 
      error.var)
    postvar0 = postbeta0$var
    postmean0 = postbeta0$mean
    postvar1 = postbeta1$var
    postmean1 = postbeta1$mean
  }
  
  # collect first new point in the stage
  D = matrix(rep(NA, N2 * dimX), N2, dimX)
  if(batch.idx == 1){
    w_candidates = apply(candidates, 1, FUN = function(x) WNlmvs(
      x, model0, model1, postmean0, postmean1, postvar0, postvar1, error.var))
    which_opt_w = which.max(w_candidates)
    x_opt_w = candidates[which_opt_w, , drop = FALSE]
    is_x_max_in_initD = any(apply(
      initD, 1, function(x) 
        isTRUE(all.equal(unname(x), unname(as.numeric(x_opt_w))))))
  } else{
    is_x_max_in_initD = TRUE
  }
  # Update set of design points (D) and plot new point
  if(is_x_max_in_initD){
    # Find minimizer of f_min
    f_min_candidates = apply(candidates, 1, function(x) f_min_vs(
      x, initD, model0, model1, postmean0, postmean1, postvar0, postvar1, 
      error.var, p, k, alpha)) 
    which_opt_f = which.min(f_min_candidates)
    x_opt_f = candidates[which_opt_f, , drop = FALSE]
    D[1, ] = x_opt_f
  } else{
    D[1, ] = x_opt_w
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find minimizer of f_min
      f_min_candidates = apply(candidates, 1, function(x) f_min_vs(
        x, rbind(initD, D[1:(i - 1), , drop = FALSE]), model0, model1, 
        postmean0, postmean1, postvar0, postvar1, error.var, p, k, alpha))
      which_opt_f = which.min(f_min_candidates)
      x_opt_f = candidates[which_opt_f, , drop = FALSE]
      # Update set of design points (D) and plot new point
      D[i, ] = x_opt_f
    }
  }
  
  return(list("initD" = initD, "addD" = D, "D" = rbind(initD, D)))
}
