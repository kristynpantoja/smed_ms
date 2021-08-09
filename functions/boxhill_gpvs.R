# using joint distribution of y | x 
Evidence_gpvs = function(y, x, model){
  null_mean_vec = rep(0, length(y))
  if(is.null(model$measurement.var)){
    K_obs = getCov(
      X1 = x[, model$indices, drop = FALSE], 
      X2 = x[, model$indices, drop = FALSE], type = model$type, l = model$l, 
      p = model$p, signal.var = model$signal.var)
  } else{
    K_obs = getCov(
      X1 = x[, model$indices, drop = FALSE], 
      X2 = x[, model$indices, drop = FALSE], type = model$type, l = model$l, 
      p = model$p, signal.var = model$signal.var) + 
      model$measurement.var * diag(length(y))
  }
  evidence = dmvnorm(
    y, mean = null_mean_vec, sigma = K_obs, log = FALSE)
  return(evidence)
}

# compute Box-Hill discrimination criterion, D, for GP case
BHDgpvs_m2 = function(
  y, # vector observed responses, at points x
  x, # vector of observed points (with data y)
  post.probs,
  candidate, # new point (potentially next point x.n, where no y.n is observed yet)
  # since these points don't have special covariates like in linear case, 
  #   they are the same for both model.i and model.j
  model.i, # type of covariance function type.i, 
  #   length-scale parameter l.i
  model.j # type of covariance function type.j, 
  #   length-scale parameter l.j
){
  # posterior predictive distributions
  candidate = t(as.matrix(candidate))
  x.i = x[, model.i$indices, drop = FALSE]
  cand.i = candidate[, model.i$indices, drop = FALSE]
  x.j = x[, model.j$indices, drop = FALSE]
  cand.j = candidate[, model.j$indices, drop = FALSE]
  pred.i = getGPPredictive(
    x = cand.i, x.input = x.i, y.input = y, type = model.i$type, l = model.i$l, 
    p = model.i$p, signal.var = model.i$signal.var, 
    measurement.var = model.i$measurement.var)
  pred.j = getGPPredictive(
    x = cand.j, x.input = x.j, y.input = y, type = model.j$type, l = model.j$l,
    p = model.j$p, signal.var = model.j$signal.var, 
    measurement.var = model.j$measurement.var)
  # evaluate criterion D
  KLij = KLN(pred.i$pred_mean, pred.i$pred_var, 
             pred.j$pred_mean, pred.j$pred_var, dim = 1)
  KLji = KLN(pred.j$pred_mean, pred.j$pred_var, 
             pred.i$pred_mean, pred.i$pred_var, dim = 1)
  bhd = prod(post.probs) * (KLij + KLji)
  return(bhd)
}

BHDgpvs_m2_testing = function(y, x, post.probs, candidate, model.i, model.j){
  # posterior predictive distributions
  x.i = x[, model.i$indices, drop = FALSE]
  cand.i = candidate[, model.i$indices, drop = FALSE]
  x.j = x[, model.j$indices, drop = FALSE]
  cand.j = candidate[, model.j$indices, drop = FALSE]
  pred.i = getGPPredictive(
    x = cand.i, x.input = x.i, y.input = y, type = model.i$type, l = model.i$l, 
    p = model.i$p, signal.var = model.i$signal.var, 
    measurement.var = model.i$measurement.var)
  pred.j = getGPPredictive(
    x = cand.j, x.input = x.j, y.input = y, type = model.j$type, l = model.j$l,
    p = model.j$p, signal.var = model.j$signal.var, 
    measurement.var = model.j$measurement.var)
  # evaluate criterion D
  KLij = KLN(pred.i$pred_mean, pred.i$pred_var, 
             pred.j$pred_mean, pred.j$pred_var, dim = 1)
  KLji = KLN(pred.j$pred_mean, pred.j$pred_var, 
             pred.i$pred_mean, pred.i$pred_var, dim = 1)
  bhd = prod(post.probs) * (KLij + KLji)
  return(list("bhd" = bhd, "KLij" = KLij, "KLji" = KLji))
}

# obtain the next n design points (fully sequential)
BHgpvs_m2 = function(
  y = NULL, # preliminary data response
  x = NULL, # preliminary data input
  x.idx = NULL, 
  prior.probs = rep(1 / 2, 2), # prior probabilities for H0 and H1
  model0, 
  model1, 
  n, # number of new points
  candidates, # domain over which the function is evaluated 
  function.values, # true function values, evaluated over the domain
  seed = NULL
  # stopping.type = "tryCatch"
){
  if(!is.null(seed)) set.seed(seed)
  
  # check preliminary data
  if(is.null(x) & !is.null(y)){ # x is null, y is not null
    stop("BHgpvs_m2 : preliminary y  is given, but not corresponding x")
  } else if(is.null(y) & !is.null(x)){ # x is not null, y is null (get y)
    y = function.values[x.idx]
  } else if(is.null(x) & is.null(y)){ # both x and y are null, then us BH method
    stop("BHgpvs_m2: need input data, at least x!")
  } else{
    if(!("matrix" %in% class(x))) stop("BHgpvs_m2: x isn't a matrix!")
    initN = dim(x)[1]
    if(length(y) != initN) stop("BHgpvs_m2: length of y does not match length of initial input data, x")
    
  }
  
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior prob with current data
    prior.probs = prior.probs, 
    evidences = c(
      Evidence_gpvs(y, x, model0),
      Evidence_gpvs(y, x, model1)
    )
  )
  # get new data
  x.new.idx = rep(NA, n)
  x.new = matrix(NA, nrow = n, ncol = ncol(candidates))
  y.new = rep(NA, n)
  post.probs.mat = matrix(NA, nrow = n + 1, ncol = 2, byrow = TRUE)
  post.probs.mat[1, ] = post.probs0
  post.probs.cur = post.probs0
  y.cur = y
  x.cur = x
  for(i in 1:n){
    bhd_seq = apply(candidates, 1, FUN = function(x) BHDgpvs_m2(
      y = y.cur, x = x.cur, post.probs = post.probs.cur, candidate = x, 
      model.i = model0, model.j = model1))
    if(!all(!is.nan(bhd_seq))){
      warning("Warning in BHgp_m2() : There were NaNs in Box & Hill criterion
              evaluation over the candidate set!!")
    }
    if(all(is.nan(bhd_seq))){
      warning("Warning in BHgp_m2() : ***ALL*** NaNs in Box & Hill criterion
              evaluated over the candidate set!!")
    }
    # get new point
    x.new.idx[i] = which.max(bhd_seq)
    x.new[i, ] = candidates[x.new.idx[i], , drop = FALSE]
    y.new[i] = function.values[x.new.idx[i]]
    # update current information
    x.cur = rbind(x.cur, x.new[i])
    y.cur = c(y.cur, y.new[i])
    post.probs.cur = getHypothesesPosteriors(
      prior.probs = post.probs.cur, 
      evidences = c(
        Evidence_gpvs(y.cur, x.cur, model0),
        Evidence_gpvs(y.cur, x.cur, model1)
      )
    )
    post.probs.mat[i + 1, ] = post.probs.cur
    
    # if(!(stopping.type == 1 & stopping.type == "tryCatch")){
    # check post.probs.cur -- if either is NaN, stop
    if(sum(is.nan(post.probs.cur)) > 0) break ## don't do for "tryCatch"
    # check post.probs.cur -- if either equals 0 or 1
    if(sum(post.probs.cur %in% c(0, 1)) > 0) break  ## don't do for "tryCatch"
    # }
  }
  return(list(
    x = x, 
    x.idx = x.idx, 
    y = y, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    post.probs = post.probs.mat, 
    function.values = function.values
  ))
}
