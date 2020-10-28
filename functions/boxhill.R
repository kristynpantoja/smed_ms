# D, criterion in Box and Hill 1967
#   for y ~ N(eta, sigma^2)
#     eta = E[y] = function that is linear (or approx linear) in parameters
#   and for m models

# get posterior distribution for beta
getPosteriorDistr = function(y, X, beta_prior_mean, beta_prior_var, var_e){
  # arguments:
  # y : response
  # X : design matrix
  
  # calculate the posterior variance (a matrix)
  beta.var.inv = NA
  # if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
  #   if(diagPrior == TRUE){
  #     beta.var.inv = (1 / beta_prior_var)
  #   } else{
  #     beta.var.inv = solve(beta_prior_var)
  #   }
  # } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
  #   if(diagPrior == TRUE){
  #     beta.var.inv = diag(1/diag(beta_prior_var))
  #   } else{
  #     beta.var.inv = solve(beta_prior_var)
  #   }
  # }
  if(diagPrior == TRUE){
    if(is.null(dim(X)) | (dim(X)[2] == 1)){
      beta.var.inv = (1 / beta_prior_var)
    } else{
      beta.var.inv = diag(1/diag(beta_prior_var))
    }
  } else{
    beta.var.inv = solve(beta_prior_var)
  }
  postvar = var_e * solve(crossprod(X) + var_e * beta.var.inv)
  
  postmean = (1 / var_e) * postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(mean = postmean, var = postvar)
}

# p_i, the probability density function of the nth observation under model i
getEvidence = function(
  y, # observed value
  sigma.t, # posterior predictive response from ith model
  yhat.i, # known sd
  sigma.i # posterior predictive variance of tilde.y.i
){
  dnorm(y, yhat.i, sigma.t + sigma.i)
}

getPosteriorProbs = function(
  prior.probs, # posterior probability Pi_{i, n - 1}
  evidences # probability density function of nth observation under model i
){
  # check that the two vectors are equal length
  if(length(prior.probs) != length(evidences)){
    stop("number of model priors doesnt match number of model evidences!")
  }
  
  # calculate vector of posterior probabilities
  m = length(prior.probs)
  unnormalized.posterior.probs = prior.probs * evidences
  posterior.probs = unnormalized.posterior.probs / 
    sum(unnormalized.posterior.probs)
  return(posterior.probs)
}

# for m = 2 models, i and j:
BHDlinear_pair_criterion = function(
  sigma.t, # data, true sd
  Pi.i, yhat.i, sigma.i, # ith model info
  Pi.j, yhat.j, sigma.j # jth model info
  ){
  # calculate pair
  inner.sigma.term = (sigma.i^2 - sigma.j^2)^2 / 
    ((sigma.t^2 + sigma.i^2) * (sigma.t^2 + sigma.j^2))
  inner.y.term = (yhat.i - yhat.j)^2 * 
    ((sigma.t^2 + sigma.i^2)^(-1) + (sigma.t^2 + sigma.j^2)^(-1))
  pair.sumterm = 0.5 * Pi.i * Pi.j * (inner.sigma.term + inner.y.term)
  return(pair.sumterm)
}

BHDlinear_pair = function(
  y, # observed responses
  sigma.t, 
  prior.probs, # 2-vector
  model.i, # design matrix of covariates for observed points X.i, 
  #   covariates of new point X.i.n, beta.prior.mu.i, beta.prior.var.i
  model.j # design matrix of covariates for observed points X.j, 
  #   covariates of new point X.j.n, beta.prior.mu.j, beta.prior.var.j
){
  # posterior predictives : yhat.i, yhat.j, sigma.i, sigma.j
  postdistr.i = getPosteriorDistr(
    y, model.i$X, model.i$beta.prior.mu, model.i$beta.prior.var, sigma.t
    )
  yhat.i = model.i$X.i.n %*% postdistr.i$mean
  sigma.i = diag(postdistr.i$var)
  postdistr.j = getPosteriorDistr(
    y, model.j$X, model.j$beta.prior.mu, model.j$beta.prior.var
    )
  yhat.j = model.j$X.j.n %*% postdistr.j$mean
  sigma.j = diag(postdistr.j$var)
  # evidences
  evidences = c(
    getEvidence(y, sigma.t, yhat.i, sigma.i), 
    getEvidence(y, sigma.t, yhat.j, sigma.j)
  )
  # posterior probs : Pi.i, Pi.j
  new.prior.probs = getPosteriorProbs(prior.probs, evidences)
  Pi.i = new.prior.probs[1]
  Pi.j = new.prior.probs[2]
  # evaluate criterion D
  BHD = BHDlinear_pair_criterion(sigma.t, Pi.i, yhat.i, sigma.i, Pi.j, yhat.j, sigma.j)
  return(BHD)
}