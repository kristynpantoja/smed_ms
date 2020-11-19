# logjointlik = function(x.new, y.new, x.old, y.old, type_arg, l_arg){
#   Sigma11 = getCov(x.old, x.old, type_arg, l_arg)
#   Sigma22 = getCov(x.new, x.new, type_arg, l_arg)
#   Sigma21 = getCov(x.old, x.new, type_arg, l_arg)
#   joint_var = rbind(cbind(Sigma11, Sigma21), cbind(t(Sigma21), Sigma22))
#   y = c(y.old, y.new)
#   return(dmvnorm(y, mean = rep(0, length(y)), sigma = joint_var, log = TRUE))
# }

# get (joint) evidence
# should I be using posterior predictive for evidence, instead? ##############################################################################
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
  for(k in 1:length(candidate)){
    KLij[k] = KLN(pred.i$pred_mean[k], pred.i$pred_var[k], 
                    pred.j$pred_mean[k], pred.j$pred_var[k], dim = 1)
    KLji[k] = KLN(pred.j$pred_mean[k], pred.j$pred_var[k], 
                    pred.i$pred_mean[k], pred.i$pred_var[k], dim = 1)
  }
  bhd = prod(post.probs) * (KLij + KLji)
  return(bhd)
}

# obtain the next n design points (fully sequential)
BHgp_m2 = function(
  y, # preliminary data response
  x, # preliminary data input
  prior.probs = rep(1 / 2, 2), # prior probabilities for H0 and H1
  model0, 
  model1, 
  n, # number of new points
  candidates, # domain over which the function is evaluated 
  function.values, # true function values, evaluated over the domain
  nugget = NULL
){
  # make sure number of elements in initD matches number of elements in y
  if(length(x) != length(y)){
    stop("Error in BHgp_m2() : length of y input doesn't match length of x 
         input!!")
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
  post.probs.mat = matrix(rep(post.probs0, n + 1), nrow = n + 1, ncol = 2, 
                          byrow = TRUE)
  post.probs.cur = post.probs0
  y.cur = y
  x.cur = x
  for(i in 1:n){
    # evaluate criterion over x_seq
    bhd_seq = BHDgp_m2(
      y.cur, x.cur, post.probs.cur, candidates, model0, model1, nugget)
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
  }
  return(list(
    x.new.idx = x.new.idx,
    x.new = x.new,
    y.new = y.new,
    post.probs = post.probs.mat
  ))
}
