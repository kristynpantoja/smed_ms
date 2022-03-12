# objective function
objective_newq_seqmed = function(
  candidate, D, postmean0, postmean1, postvar0, postvar1, model0, model1, 
  error.var, p = 1, k = NULL, alpha = 1){
  if(is.null(k)) k = 4 * p
  # result0 = q_seqmed(
  #   candidate, postmean0, postmean1, postvar0, postvar1, model0, model1,
  #   error.var, p, alpha)^k *
  #   sum(sapply(
  #     D,
  #     function(x_i) (
  #       q_seqmed(x_i, postmean0, postmean1, postvar0, postvar1,
  #                model0, model1, error.var, p, alpha) /
  #         sqrt((x_i - candidate)^2)
  #     )^k))
  
  q_cand = q_seqmed(
    x = candidate, postmean0 = postmean0, postmean1 = postmean1, 
    postvar0 = postvar0, postvar1 = postvar1, model0 = model0, model1 = model1, 
    error.var = error.var, p = p, alpha = alpha)
  q_D = sapply(
    D, 
    function(x_i)
      q_seqmed(
        x = x_i, postmean0 = postmean0, postmean1 = postmean1, 
        postvar0 = postvar0, postvar1 = postvar1, 
        model0 = model0, model1 = model1, error.var = error.var, p = p, 
        alpha = alpha) / 
      sqrt((x_i - candidate)^2)
  )
  result = q_cand^k * sum(q_D^k)
  return(result)
}

# batch of seqN new points
SeqMED_newq_batch = function(
  initD, y, model0, model1, error.var, 
  N2 = 1, numCandidates = 10^5, k = NULL, 
  xmin = -1, xmax = 1, p = 1, alpha = 1, genCandidates = 1, candidates = NULL, 
  keep_trying_alpha = keep_trying_alpha, 
  prints = FALSE, save_objectives = FALSE, batch.idx = 1
){
  if(is.null(k)) k = 4 * p
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # posterior distributions of beta
  postbeta0 = getBetaPosterior(
    y = y, X = model0$designMat(initD), model0$beta.mean, model0$beta.var, 
    error.var)
  postbeta1 = getBetaPosterior(
    y = y, X = model1$designMat(initD), model1$beta.mean, model1$beta.var, 
    error.var)
  postvar0 = postbeta0$var
  postmean0 = postbeta0$mean
  postvar1 = postbeta1$var
  postmean1 = postbeta1$mean
  
  # check if any points in initD give Wasserstein distance of 0
  #   if there are any such points, remove them so that TPE does not explode
  #   (since 1/0 in q)
  w_initD = sapply(initD, FUN = function(x) WNlm(
    x, postmean0, postmean1, postvar0, postvar1, model0, model1, error.var))
  if(length(which(w_initD == 0)) != 0){
    initD = initD[-which(w_initD == 0)]
    y = y[-which(w_initD == 0)]
    w_initD = w_initD[-which(w_initD == 0)]
    # recalculate posterior distributions of beta
    postbeta0 = getBetaPosterior(
      y = y, X = model0$designMat(initD), model0$beta.mean, model0$beta.var, 
      error.var)
    postbeta1 = getBetaPosterior(
      y = y, X = model1$designMat(initD), model1$beta.mean, model1$beta.var, 
      error.var)
    postvar0 = postbeta0$var
    postmean0 = postbeta0$mean
    postvar1 = postbeta1$var
    postmean1 = postbeta1$mean
  }
  
  # generate candidate points, if necessary
  if(is.null(candidates)){
    if(genCandidates == 1){
      candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    } else if(genCandidates == 2){
      candidates = sort(runif(numCandidates, min = xmin, max = xmax))
    }
  }
  
  # collect first new point in the stage
  D = rep(NA, N2)
  if(batch.idx == 1){
    w_candidates = sapply(candidates, function(x) WNlm(
      x, postmean0, postmean1, postvar0, postvar1, model0, model1, error.var))
    which_opt_w = which.max(w_candidates)
    x_opt_w = candidates[which_opt_w]
    is_x_max_in_initD = any(sapply(initD, function(x) x == x_opt_w))
    if(save_objectives){
      objectives = w_candidates
    }
  } else{
    is_x_max_in_initD = TRUE
  }
  # Update set of design points (D) and plot new point
  if(is_x_max_in_initD){
    # Find minimizer of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) objective_newq_seqmed(
        candidate = x, D = initD, 
        postmean0 = postmean0, postmean1 = postmean1, 
        postvar0 = postvar0, postvar1 = postvar1, 
        model0 = model0, model1 = model1, error.var = error.var, 
        p = p, k = k, alpha = alpha))
    if(prints){
      if(any(f_min_candidates == Inf)){
        print(paste0(
          "SeqMED_newq_batch: objective function evaluates to Inf for ",
          as.character(sum(f_min_candidates == Inf)), 
          " candidates"
        ))
      }
    }
    if(all(f_min_candidates == Inf)){
      warning("SeqMED_newq_batch: objective function evaluates to Inf for all candidates")
      if(keep_trying_alpha){
        while(alpha > 0 & all(f_min_candidates == Inf)){
          alpha = alpha - 1
          f_min_candidates = sapply(
            candidates, 
            function(x) objective_newq_seqmed(
              candidate = x, D = initD, 
              postmean0 = postmean0, postmean1 = postmean1, 
              postvar0 = postvar0, postvar1 = postvar1, 
              model0 = model0, model1 = model1, error.var = error.var, 
              p = p, k = k, alpha = alpha))
        }
        if(alpha <= 0){
          stop("SeqMED_newq_batch: All alpha > 0 result in the objective function evaluating to Inf")
        }
      }
    }
    if(save_objectives){
      objectives = f_min_candidates
    }
    which_opt_f = which.min(f_min_candidates)
    x_opt_f = candidates[which_opt_f]
    D[1] = x_opt_f
  } else{
    D[1] = x_opt_w
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find minimizer of f_min
      f_min_candidates = sapply(
        candidates, 
        function(x) objective_newq_seqmed(
          candidate = x, D = c(initD, D[1:(i - 1)]), 
          postmean0 = postmean0, postmean1 = postmean1, 
          postvar0 = postvar0, postvar1 = postvar1, 
          model0 = model0, model1 = model1, error.var = error.var, 
          p = p, k = k, alpha = alpha))
      if(prints){
        if(any(f_min_candidates == Inf)){
          print(paste0(
            "SeqMED_newq_batch: objective function evaluates to Inf for ",
            as.character(sum(f_min_candidates == Inf)), 
            " candidates"
          ))
        }
      }
      if(all(f_min_candidates == Inf)){
        warning("SeqMED_newq_batch: objective function evaluates to Inf for all candidates")
        if(keep_trying_alpha){
          while(alpha > 0  & all(f_min_candidates == Inf)){
            alpha = alpha - 1
            f_min_candidates = sapply(
              candidates, 
              function(x) objective_newq_seqmed(
                candidate = x, D = c(initD, D[1:(i - 1)]), 
                postmean0 = postmean0, postmean1 = postmean1, 
                postvar0 = postvar0, postvar1 = postvar1, 
                model0 = model0, model1 = model1, error.var = error.var, 
                p = p, k = k, alpha = alpha))
          }
        }
      }
      which_opt_f = which.min(f_min_candidates)
      x_opt_f = candidates[which_opt_f]
      # Update set of design points (D) and plot new point
      D[i] = x_opt_f
    }
  }
  if(prints){
    print(paste0("alpha=", as.character(alpha)))
  }
  result = list(
    initD = initD, addD = D, updatedD = c(initD, D), alpha = alpha
  )
  if(save_objectives){
    result$objectives = objectives
  }
  return(result)
}
