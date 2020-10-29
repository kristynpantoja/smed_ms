# logjointlik = function(x.new, y.new, x.old, y.old, type_arg, l_arg){
#   Sigma11 = getCov(x.old, x.old, type_arg, l_arg)
#   Sigma22 = getCov(x.new, x.new, type_arg, l_arg)
#   Sigma21 = getCov(x.old, x.new, type_arg, l_arg)
#   joint_var = rbind(cbind(Sigma11, Sigma21), cbind(t(Sigma21), Sigma22))
#   y = c(y.old, y.new)
#   return(dmvnorm(y, mean = rep(0, length(y)), sigma = joint_var, log = TRUE))
# }

getEvidenceGP = function(
  y, x, x.n,
  model
){
  # get (joint) evidences
  evidence = dmvnorm(
    y, mean = rep(0, length(y)), 
    sigma = getCov(c(x, x.n), model$type, model$l), log = FALSE)
  return(evidence)
}

getKLMVN(
  mu1, Sigma1, 
  mu2, Sigma2
){
  if(length(mu1) != length(mu2)) stop("Error in getKLMVN : mu1 and mu2 are not same length!")
  if(length(mu1) != dim(Sigma1)) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  if(length(mu2) != dim(Sigma2)) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  d = length(mu1)
  log.det.term = log(det(Sigma2)) - log(det(Sigma1))
  tr.term = sum(diag(solve(Sigma2, Sigma1))) * solve(Sigma2, mu2 - mu1)
  0.5 * (log.det.term - d + tr.term + quadratic.term)
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
  pred.i = getPredDistrSeq(x.n, x, y, model.i$type, model.i$l, nugget)
  pred.j = getPredDistrSeq(x.n, x, y, model.j$type, model.j$l, nugget)
  yhat.i = pred.i$pred_mean
  yhat.j = pred.j$pred_mean
  # evidences
  evidences = c(
    getEvidenceGP(y, x, x.n, model.i), 
    getEvidenceGP(y, x, x.n, model.j)
  )
  # posterior probs : Pi.i, Pi.j
  new.prior.probs = getPosteriorProbs(prior.probs, evidences)
  Pi.i = new.prior.probs[1]
  Pi.j = new.prior.probs[2]
  # evaluate criterion D
  KLij = getKLGP(pred.i$pred_mean, pred.i$pred_var, 
                 pred.j$pred_mean, pred.j$pred_var)
  KLji = getKLGP(pred.j$pred_mean, pred.j$pred_var, 
                 pred.i$pred_mean, pred.i$pred_var)
  BHD = Pi.i * Pi.j * KLij + KLji
  return(BHD)
}