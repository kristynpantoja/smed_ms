# using joint distribution of y | x 
#   -- should I be using the posterior predictive, instead?
Evidence_gp = function(y, x, model){
  evidence = dmvnorm(
    y, mean = rep(0, length(y)), 
    sigma = getCov(x, x, model$type, model$l), log = FALSE)
  return(evidence)
}

# compute Box-Hill discrimination criterion, D, for GP case
BHDgp_m2 = function(
  y, # vector observed responses, at points x
  x, # vector of observed points (with data y)
  post.probs,
  candidate, # new point (potentially next point x.n, where no y.n is observed yet)
  # since these points don't have special covariates like in linear case, 
  #   they are the same for both model.i and model.j
  model.i, # type of covariance function type.i, 
  #   length-scale parameter l.i
  model.j, # type of covariance function type.j, 
  #   length-scale parameter l.j
  nugget = NULL
){
  # posterior predictive distributions
  pred.i = getGPPredictive(candidate, x, y, model.i$type, model.i$l, nugget)
  pred.j = getGPPredictive(candidate, x, y, model.j$type, model.j$l, nugget)
  # evaluate criterion D
  KLij = rep(NA, length(candidate))
  KLji = rep(NA, length(candidate))
  KLij = KLN(pred.i$pred_mean, pred.i$pred_var, 
             pred.j$pred_mean, pred.j$pred_var, dim = 1)
  KLji = KLN(pred.j$pred_mean, pred.j$pred_var, 
             pred.i$pred_mean, pred.i$pred_var, dim = 1)
  bhd = prod(post.probs) * (KLij + KLji)
  return(bhd)
}

# obtain the next n design points (fully sequential)
BHgp_m2 = function(
  y = NULL, # preliminary data response
  x = NULL, # preliminary data input
  x.idx = NULL, 
  prior.probs = rep(1 / 2, 2), # prior probabilities for H0 and H1
  model0, 
  model1, 
  n, # number of new points
  candidates, # domain over which the function is evaluated 
  function.values, # true function values, evaluated over the domain
  nugget = NULL,
  seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  
  # check preliminary data
  if(is.null(x) & !is.null(y)){ # x is null, y is not null
    stop("BH_m2 : preliminary y  is given, but not corresponding x")
  } else if(is.null(y) & !is.null(x)){ # x is not null, y is null (get y)
    y = function.values[x.idx]
  } else if(is.null(x) & is.null(y)){ # both x and y are null, then us BH method
    stop("BHgp_m2: need input data, at least x!")
  } else{
    if(length(x) != length(y)){
      stop("BH_m2 : length of preliminary x and y don't match!")
    }
  }
  
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior prob with current data
    prior.probs = prior.probs, 
    evidences = c(
      Evidence_gp(y, x, model0),
      Evidence_gp(y, x, model1)
    )
  )
  # get new data
  x.new.idx = rep(NA, n)
  x.new = rep(NA, n)
  y.new = rep(NA, n)
  post.probs.mat = matrix(NA, nrow = n + 1, ncol = 2, byrow = TRUE)
  post.probs.mat[1, ] = post.probs0
  post.probs.cur = post.probs0
  y.cur = y
  x.cur = x
  for(i in 1:n){
    # evaluate criterion over x_seq
    bhd_seq = sapply(candidates, FUN = function(x) BHDgp_m2(
      y.cur, x.cur, post.probs.cur, x, model0, model1, nugget))
    if(!all(!is.nan(bhd_seq))){
      warning("Warning in BHgp_m2() : There were NaNs in Box & Hill criterion 
              evaluation over the candidate set!!")
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
    # check post.probs.cur
    if(all(post.probs.cur %in% c(0, 1))) break
  }
  return(list(
    x.new.idx = x.new.idx,
    x.new = x.new,
    y.new = y.new,
    post.probs = post.probs.mat
  ))
}
