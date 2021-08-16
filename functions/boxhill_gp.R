# using joint distribution of y | x 
#   -- should I be using the posterior predictive, instead?
Evidence_gp = function(y, x, model){
  null_mean_vec = rep(0, length(y))
  if(is.null(model$measurement.var)){
    K_obs = getCov(
      X1 = x, X2 = x, type = model$type, l = model$l, p = model$p, 
      signal.var = model$signal.var)
  } else{
    K_obs = getCov(
      X1 = x, X2 = x, type = model$type, l = model$l, p = model$p, 
      signal.var = model$signal.var) + model$measurement.var * diag(length(y))
  }
  evidence = dmvnorm(
    y, mean = null_mean_vec, sigma = K_obs, log = FALSE)
  return(evidence)
}

# compute Box-Hill discrimination criterion, D, for GP case
BHDgp_m2 = function(
  y.in, # vector observed responses, at points x
  x.in, # vector of observed points (with data y)
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
  pred.i = getGPPredictive(
    x = candidate, x.input = x.in, y.input = y.in, type = model.i$type, l = model.i$l, 
    p = model.i$p, signal.var = model.i$signal.var, 
    measurement.var = model.i$measurement.var)
  pred.j = getGPPredictive(
    x = candidate, x.input = x.in, y.input = y.in, type = model.j$type, l = model.j$l,
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

# BHDgp_m2_testing = function(
#   y.in, x.in, post.probs, candidate, model.i, model.j){
#   # posterior predictive distributions
#   pred.i = getGPPredictive(
#     x = candidate, x.input = x.in, y.input = y.in, type = model.i$type, 
#     l = model.i$l, p = model.i$p, signal.var = model.i$signal.var,
#     measurement.var = model.i$measurement.var)
#   pred.j = getGPPredictive(
#     x = candidate, x.input = x.in, y.input = y.in, type = model.j$type, 
#     l = model.j$l, p = model.j$p, signal.var = model.j$signal.var,
#     measurement.var = model.j$measurement.var)
#   # evaluate criterion D
#   KLij = KLN(pred.i$pred_mean, pred.i$pred_var, 
#              pred.j$pred_mean, pred.j$pred_var, dim = 1)
#   KLji = KLN(pred.j$pred_mean, pred.j$pred_var, 
#              pred.i$pred_mean, pred.i$pred_var, dim = 1)
#   bhd = prod(post.probs) * (KLij + KLji)
#   return(list("bhd" = bhd, "KLij" = KLij, "KLji" = KLji))
# }

# obtain the next n design points (fully sequential)
BHgp_m2 = function(
  y.in = NULL, # preliminary data response
  x.in = NULL, # preliminary data input
  x.in.idx = NULL, 
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
  if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
    stop("BH_m2 : preliminary y.in is given, but not corresponding x.in")
  } else if(is.null(y.in) & !is.null(x.in)){ # x.in is not null, y.in is null (get y.in)
    y.in = function.values[x.in.idx]
  } else if(is.null(x.in) & is.null(y.in)){ # both x and y are null, then us BH method
    stop("BHgp_m2: need input data, at least x.in!")
  } else{
    if(length(x.in) != length(y.in)){
      stop("BH_m2 : length of preliminary x.in and y.in don't match!")
    }
  }
  
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior prob with current data
    prior.probs = prior.probs, 
    evidences = c(
      Evidence_gp(y.in, x.in, model0),
      Evidence_gp(y.in, x.in, model1)
    )
  )
  # get new data
  x.new.idx = rep(NA, n)
  x.new = rep(NA, n)
  y.new = rep(NA, n)
  post.probs.mat = matrix(NA, nrow = n + 1, ncol = 2, byrow = TRUE)
  post.probs.mat[1, ] = post.probs0
  post.probs.cur = post.probs0
  y.cur = y.in
  x.cur = x.in
  for(i in 1:n){
    bhd_seq = sapply(candidates, FUN = function(x) BHDgp_m2(
      y.in = y.cur, x.in = x.cur, post.probs = post.probs.cur, candidate = x, 
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
    x.new[i] = candidates[x.new.idx[i]]
    y.new[i] = function.values[x.new.idx[i]]
    # update current information
    x.cur = c(x.cur, x.new[i])
    y.cur = c(y.cur, y.new[i])
    post.probs.cur = getHypothesesPosteriors(
      prior.probs = post.probs.cur, 
      evidences = c(
        Evidence_gp(y.cur, x.cur, model0),
        Evidence_gp(y.cur, x.cur, model1)
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
    x.in = x.in, 
    x.in.idx = x.in.idx, 
    y.in = y.in, 
    x.new = x.new,
    x.new.idx = x.new.idx,
    y.new = y.new,
    post.probs = post.probs.mat, 
    function.values = function.values
  ))
}
