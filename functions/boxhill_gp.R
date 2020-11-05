# logjointlik = function(x.new, y.new, x.old, y.old, type_arg, l_arg){
#   Sigma11 = getCov(x.old, x.old, type_arg, l_arg)
#   Sigma22 = getCov(x.new, x.new, type_arg, l_arg)
#   Sigma21 = getCov(x.old, x.new, type_arg, l_arg)
#   joint_var = rbind(cbind(Sigma11, Sigma21), cbind(t(Sigma21), Sigma22))
#   y = c(y.old, y.new)
#   return(dmvnorm(y, mean = rep(0, length(y)), sigma = joint_var, log = TRUE))
# }

getEvidenceGP = function(
  y, x,
  model
){
  # get (joint) evidences
  evidence = dmvnorm(
    y, mean = rep(0, length(y)), 
    sigma = getCov(x, x, model$type, model$l), log = FALSE)
  return(evidence)
}

getKLMVN = function(
  mu1, Sigma1, 
  mu2, Sigma2
){
  if(length(mu1) != length(mu2)) stop("Error in getKLMVN : mu1 and mu2 are not same length!")
  if(length(mu1) != dim(Sigma1)[1]) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  if(length(mu2) != dim(Sigma2)[1]) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  d = length(mu1)
  log.det.term = log(det(Sigma2)) - log(det(Sigma1))
  trace.term = sum(diag(solve(Sigma2, Sigma1))) 
  quadratic.term = t(mu2 - mu1) %*% solve(Sigma2, mu2 - mu1)
  KL = 0.5 * (log.det.term - d + trace.term + quadratic.term)
  return(KL)
}

getKLN1 = function(
  mu1, sigma1, 
  mu2, sigma2
){
  KL = log(sigma2) - log(sigma1) + ((sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2)) - 0.5
}

BHDgp_m2 = function(
  y, # observed response, at points x
  post.probs,
  x, # vector of observed points (with data y)
  x.n, # new point (no y.n observed yet)
  # since these points don't have special covariates like in linear case, 
  #   they are the same for both model.i and model.j
  model.i, # type of covariance function type.i, 
  #   length-scale parameter l.i
  model.j, # type of covariance function type.j, 
  #   length-scale parameter l.j
  nugget = NULL
){
  # posterior predictive distributions
  pred.i = getGPPredictive(x.n, x, y, model.i$type, model.i$l, nugget)
  pred.j = getGPPredictive(x.n, x, y, model.j$type, model.j$l, nugget)
  yhat.i = pred.i$pred_mean
  yhat.j = pred.j$pred_mean
  # evidences
  evidences = c(
    getEvidenceGP(y, x, model.i), 
    getEvidenceGP(y, x, model.j)
  )
  # posterior probs: Pi_{ni}, Pi_{nj}
  #   (called "prior probability," or, "probability of current data")
  # post.probs = getPosteriorProbs(prior.probs, evidences)
  # evaluate criterion D
  KLij = rep(NA, length(x.n))
  KLji = rep(NA, length(x.n))
  for(k in 1:length(x.n)){
    KLij[k] = getKLN1(pred.i$pred_mean[k], diag(pred.i$pred_var)[k], 
                    pred.j$pred_mean[k], diag(pred.j$pred_var)[k])
    KLji[k] = getKLN1(pred.j$pred_mean[k], diag(pred.j$pred_var)[k], 
                    pred.i$pred_mean[k], diag(pred.i$pred_var)[k])
  }
  BHD = prod(post.probs) * (KLij + KLji)
  return(BHD)
}

add1_BHgp_m2 = function(
  y, 
  post.probs,
  x, 
  model0, # H0 type of covariance function, length-scale parameter l
  model1, # H1 type of covariance function, length-scale parameter l
  domain, 
  nugget = NULL
  ){
  
  # make sure number of elements in initD matches number of elements in y
  if(length(y) != length(y)) stop("Error in add1_BHgp_m2() : length of y does not match length of initial input data, initD!!")
  
  # calculate posterior probabilities, given current data
  # posterior probs: Pi_{ni}, Pi_{nj}
  #   (called "prior probability," or, "probability of current data")
  # post.probs = getPosteriorProbs(prior.probs, evidences)
  
  # evaluate criterion over x_seq
  BHD_seq = BHDgp_m2(y, post.probs, x, domain, model0, model1, nugget)
  if(!all(!is.nan(BHD_seq))) warning("Warning in BHDgp_m2() : There were NaNs in Box & Hill criterion evaluation over the candidate set!!")
  # which(is.nan(BHD_seq))
  
  # get new point
  x.new.idx = which.max(BHD_seq)
  x.new = x_seq[x.new.idx]
  
  return(list(
    x.new = x.new,
    x.new.idx = x.new.idx
    ))
}

getBHgp_m2 = function(
  y, # preliminary data response
  prior.probs = rep(1 / 2, 2), # prior probabilities for H0 and H1
  x, # preliminary data
  model0, 
  model1, 
  domain, 
  function.values, # true function values, evaluated over the domain
  n, # number of new points
  nugget = NULL
){
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getPosteriorProbs( # posterior probability with current data
    prior.probs = prior.probs, 
    evidences = c(
      getEvidenceGP(y, x, model0),
      getEvidenceGP(y, x, model1)
    )
  )
  
  # get new data
  x.new.idx = rep(NA, n)
  x.new = rep(NA, n)
  y.new = rep(NA, n)
  post.probs.mat = matrix(rep(post.probs0, n + 1), nrow = n + 1, ncol = 2, byrow = TRUE)
  post.probs.cur = post.probs0
  y.cur = y
  x.cur = x
  for(i in 1:n){
    BHrun = add1_BHgp_m2(y.cur, post.probs.cur, x.cur, model0, model1, domain, nugget)
    x.new.idx[i] = BHrun$x.new.idx
    x.new[i] = BHrun$x.new
    y.new[i] = function.values[x.new.idx[i]]
    # update current information
    x.cur = c(x.cur, x.new[i])
    y.cur = c(y.cur, y.new[i])
    post.probs.cur = getPosteriorProbs(
      prior.probs = post.probs.cur, 
      evidences = c(
        getEvidenceGP(y.cur, x.cur, model0),
        getEvidenceGP(y.cur, x.cur, model1)
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
