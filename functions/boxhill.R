# D, criterion in Box and Hill 1967
#   for y ~ N(eta, sigma^2)
#     eta = E[y] = function that is linear (or approx linear) in parameters
#   and for m models

# get posterior distribution for beta
getBetaPosterior = function(
  y, X, beta_prior_mean, beta_prior_var, error.var, diagPrior = TRUE){
  # y : response
  # X : design matrix
  
  # calculate the posterior variance (a matrix)
  # beta.var.inv = solve(beta_prior_var)
  beta.var.inv = NA
  if(diagPrior == TRUE){
    if(is.null(dim(X)) | (dim(X)[2] == 1)){
      beta.var.inv = (1 / beta_prior_var)
    } else{
      beta.var.inv = diag(1 / diag(beta_prior_var))
    }
  } else{
    beta.var.inv = solve(beta_prior_var)
  }
  postvar = error.var * solve(crossprod(X) + error.var * beta.var.inv)
  postmean = (1 / error.var) * postvar %*% (t(X) %*% y + error.var
                                            * solve(beta_prior_var, beta_prior_mean)
                                            )
  return(list(mean = postmean, var = postvar))
}

# univariate
getLMPredictive = function(
  X.n, y, X, beta_prior_mean, beta_prior_var, error.var, diagPrior = TRUE){
  postbeta = getBetaPosterior(
    y, X, beta_prior_mean, beta_prior_var, error.var, diagPrior = TRUE)
  postpred = list( # posterior predictive
    mean = X.n %*% postbeta$mean, 
    var = X.n %*% postbeta$var %*% t(X.n) + error.var)
  return(postpred)
}

getHypothesesPosteriors = function(
  prior.probs, # probabilities Pi_{i, n - 1}
  evidences # probability density function of nth observation, under each model
){
  # check that the two vectors are equal length
  if(length(prior.probs) != length(evidences)){
    stop("number of model priors doesnt match number of model evidences!")
  }
  # calculate vector of posterior probabilities
  unnormalized.posterior.probs = prior.probs * evidences
  posterior.probs = unnormalized.posterior.probs / 
    sum(unnormalized.posterior.probs)
  return(posterior.probs)
}

# get evidence
Evidence_lm = function(y, x, model, error.var){
  # posterior predictive distribution
  # postpred = getLMPredictive(
  #   model$designMat(x), y, model$designMat(x), model$beta.mean, model$beta.var, error.var)
  # evidences: p_i(nth observation), where p_i is the density under model i
  # evidence = dmvnorm(y, postpred$mean, postpred$var)
  ### try marginal #################################################################
  n = length(y)
  # get mean and variance of marginal density of y, which is n-dim multivariate normal pdf
  X = model$designMat(x)
  marginaly_mean = X %*% model$beta.mean
  if(dim(X)[1] > 1){ # if X is a matrix of inputs
    marginaly_var = diag(rep(error.var, n)) + (X %*% model$beta.var %*% t(X))
  } else{ # if X is a vector for one input
    warning("X is a vector, not a matrix - is that what you expected?")
    marginaly_var = error.var + (X %*% model$beta.var %*% t(X))
  }
  evidence = dmvnorm(y, mean = marginaly_mean, sigma = marginaly_var)
  return(evidence)
}

# for m = 2 models, i and j:
# for single candidate
BHD_m2 = function(
  y, # vector observed responses, at points x
  x, # vector of observed points (with data y)
  post.probs,
  candidate, # new point (potentially next point x.n, where no y.n is observed yet)
  model.i, # DesignMat, beta.mean, beta.var -- for model i
  model.j, # DesignMat, beta.mean, beta.var -- for model j
  error.var # known error variance
){
  # posterior predictive distributions
  postpred.i = getLMPredictive(
    model.i$designMat(candidate), y, model.i$designMat(x), model.i$beta.mean, model.i$beta.var, error.var)
  postpred.j = getLMPredictive(
    model.j$designMat(candidate), y, model.j$designMat(x), model.j$beta.mean, model.j$beta.var, error.var)
  # evaluate criterion D
  ########################################################################################## treating y's indep.....
  first.term = (postpred.i$var - postpred.j$var)^2 / 
    (postpred.i$var * postpred.j$var)
  second.term = (postpred.i$mean - postpred.j$mean)^2 *
    (postpred.i$var^(-1) + postpred.j$var^(-1))
  bhd = 0.5 * prod(post.probs) * (first.term + second.term)
  return(bhd)
}

BH_m2 = function(
  y, # preliminary data response
  x, # preliminary data input
  prior.probs = rep(1 / 2, 2), # prior probabilities for H0 and H1
  model0, # DesignMat, beta.mean, beta.var -- for model0
  model1, # DesignMat, beta.mean, beta.var -- for model1
  n, # number of new points
  candidates, # set of candidate points in the domain
  true.function, # true function, linear in parameters
  error.var, # known error variance
  seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  # make sure number of elements in initD matches number of elements in y
  if(length(y) != length(x)){
    stop("Error in BH_m2() : length of y input doesn't match length of x 
         input!!")
  }
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior probability with current data
    prior.probs = prior.probs, 
    evidences = c(
      Evidence_lm(y, x, model0, error.var), 
      Evidence_lm(y, x, model1, error.var)
    )
  )
  # get new data
  x.new = rep(NA, n)
  y.new = rep(NA, n)
  post.probs.mat = matrix(rep(post.probs0, n + 1), nrow = n + 1, ncol = 2, 
                          byrow = TRUE)
  post.probs.cur = post.probs0
  y.cur = y
  x.cur = x
  for(i in 1:n){
    # evaluate criterion over x_seq
    bhd_seq = sapply(
      candidates, 
      FUN = function(x) BHD_m2(y.cur, x.cur, post.probs.cur, x, model0, model1, 
                               error.var))
    # bhd_seq = BHD_m2(
    #   y.cur, x.cur, post.probs.cur, candidates, model0, model1, error.var)
    if(!all(!is.nan(bhd_seq))){
      warning("Warning in BHDgp_m2() : There were NaNs in Box & Hill criterion 
              evaluation over the candidate set!!")
    }
    # get new point
    x.new.idx = which.max(bhd_seq)
    x.new[i] = candidates[x.new.idx]
    # y.new[i] = true.function(x.new[i])
    y.new[i] = simulateY_fromfunction(x.new[i], true.function, error.var, 
                                      numSims = 1)
    # update current information
    x.cur = c(x.cur, x.new[i])
    y.cur = c(y.cur, y.new[i])
    post.probs.cur = getHypothesesPosteriors(
      prior.probs = post.probs.cur, 
      evidences = c(
        Evidence_lm(y.cur, x.cur, model0, error.var),
        Evidence_lm(y.cur, x.cur, model1, error.var)
      )
    )
    post.probs.mat[i + 1, ] = post.probs.cur
  }
  return(list(
    x.new = x.new,
    y.new = y.new,
    post.probs = post.probs.mat
  ))
  return(BHD)
}
