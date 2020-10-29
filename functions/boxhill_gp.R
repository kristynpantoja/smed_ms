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

BHDgp_pair = function(
  y, # observed response, at points x
  prior.probs,
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
  # posterior probs : Pi.i, Pi.j
  new.prior.probs = getPosteriorProbs(prior.probs, evidences)
  Pi.i = new.prior.probs[1]
  Pi.j = new.prior.probs[2]
  # evaluate criterion D
  KLij = rep(NA, length(x.n))
  KLji = rep(NA, length(x.n))
  for(k in 1:length(x.n)){
    KLij[k] = getKLN1(pred.i$pred_mean[k], diag(pred.i$pred_var)[k], 
                    pred.j$pred_mean[k], diag(pred.j$pred_var)[k])
    KLji[k] = getKLN1(pred.j$pred_mean[k], diag(pred.j$pred_var)[k], 
                    pred.i$pred_mean[k], diag(pred.i$pred_var)[k])
  }
  BHD = Pi.i * Pi.j * (KLij + KLji)
  return(BHD)
}