# D, criterion in Box and Hill 1967
#   for y ~ N(eta, sigma^2)
#     eta = E[y] = function that is linear (or approx linear) in parameters
#   and for m models

# get posterior distribution for beta
getBetaPosterior = function(
  y, 
  X, 
  beta.mean, 
  beta.var, 
  error.var, 
  diagPrior = TRUE
){
  # y : response
  # X : design matrix
  
  # calculate the posterior variance (a matrix)
  # beta.var.inv = solve(beta.var)
  beta.var.inv = NA
  if(diagPrior == TRUE){
    if(!("matrix" %in% class(X)) | (ncol(X) == 1)){ # univariate
      beta.var.inv = (1 / beta.var)
    } else{ # multivariate case, 
      #         e.g. lm with more than 1 parameter, or lmvs with dim > 1
      beta.var.inv = diag(1 / diag(beta.var))
    }
  } else{
    beta.var.inv = solve(beta.var)
  }
  postvar = error.var * solve(crossprod(X) + error.var * beta.var.inv)
  postmean = (1 / error.var) * 
    postvar %*% (t(X) %*% y + error.var * solve(beta.var, beta.mean))
  return(list(
    mean = postmean, 
    var = postvar))
}

# univariate
getLMPredictive = function(
  X.n, 
  y, 
  X, 
  beta.mean, 
  beta.var, 
  error.var, 
  diagPrior = TRUE
){
  postbeta = getBetaPosterior(
    y, X, beta.mean, beta.var, error.var, diagPrior = TRUE)
  predmean = X.n %*% postbeta$mean
  predvar = X.n %*% postbeta$var %*% t(X.n) + error.var
  return(list( # posterior predictive
    mean = predmean, 
    var = predvar
  ))
}

getLMMarginal = function(
  X, 
  beta.mean, 
  beta.var, 
  error.var
){
  marginaly_mean = X %*% beta.mean
  if(dim(X)[1] > 1){ # if X is a matrix of inputs
    n = dim(X)[1]
    marginaly_var = diag(rep(error.var, n)) + (X %*% beta.var %*% t(X))
  } else{ # if X is a vector for one input
    marginaly_var = error.var + (X %*% beta.var %*% t(X))
  }
  return(list(
    mean = marginaly_mean, 
    var = marginaly_var
  ))
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
Evidence_lm = function(
  y, 
  x, 
  model, 
  error.var
){
  # get mean and variance of marginal density of y, which is n-dim multivariate normal pdf
  marg = getLMMarginal(model$designMat(x), model$beta.mean, model$beta.var, 
                       error.var)
  evidence = dmvnorm(y, mean = marg$mean, sigma = marg$var)
  return(evidence)
}

# for m = 2 models, i and j:
# for single candidate
BHD_m2 = function(
  y.in = NULL, # vector observed responses, at points x
  x.in = NULL, # vector of observed points (with data y)
  probs,
  candidate, # new point (potentially next point x.n, where no y.n is observed yet)
  model.i, # DesignMat, beta.mean, beta.var -- for model i
  model.j, # DesignMat, beta.mean, beta.var -- for model j
  error.var # known error variance
){
  if(is.null(x.in) & is.null(y.in)){ # use marginals
    distr.i = getLMMarginal(model.i$designMat(candidate), model.i$beta.mean, 
                            model.i$beta.var, error.var)
    distr.j = getLMMarginal(model.j$designMat(candidate), model.j$beta.mean, 
                            model.j$beta.var, error.var)
  } else{ # use posterior predictives
    distr.i = getLMPredictive(
      model.i$designMat(candidate), y.in, model.i$designMat(x.in), model.i$beta.mean, 
      model.i$beta.var, error.var)
    distr.j = getLMPredictive(
      model.j$designMat(candidate), y.in, model.j$designMat(x.in), model.j$beta.mean, 
      model.j$beta.var, error.var)
  }
  # evaluate criterion D
  first.term = (distr.i$var - distr.j$var)^2 / 
    (distr.i$var * distr.j$var)
  second.term = (distr.i$mean - distr.j$mean)^2 *
    (distr.i$var^(-1) + distr.j$var^(-1))
  bhd = 0.5 * prod(probs) * (first.term + second.term)
  return(bhd)
}

BH_m2 = function(
  y.in = NULL, # preliminary data response
  x.in = NULL, # preliminary data input
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
  
  # if no preliminary data, get the first point
  if(is.null(x.in) & !is.null(y.in)){ # x.in is null, y.in is not null
    stop("BH_m2 : preliminary y.in is given, but not corresponding x.in")
  } else if(is.null(y.in) & !is.null(x.in)){ # x.in is not null, y.in is null (generate y.in)
    y.in = simulateY_fromfunction(
      x = x.in, true.function = true.function, error.var = error.var)
  } else if(is.null(x.in) & is.null(y.in)){ # both x.in and y.in are null, then us BH method
    # evaluate criterion over x_seq
    bhd_seq = sapply(
      candidates, 
      FUN = function(x) BHD_m2(NULL, NULL, prior.probs, x, model0, model1, 
                               error.var))
    if(!all(!is.nan(bhd_seq))){
      warning("Warning in BHDgp_m2() : There were NaNs in Box & Hill criterion 
              evaluation over the candidate set!!")
    }
    # get x.in
    x.in.idx = which.max(bhd_seq)
    x.in = candidates[x.in.idx]
    # get y.in
    y.in = simulateY_fromfunction(
      x = x.in, true.function = true.function, error.var = error.var)
    n = n - 1
  } else{
    if(length(x.in) != length(y.in)){
      stop("BH_m2 : length of preliminary x.in and y.in don't match!")
    }
  }
  
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior probability with current data
    prior.probs = prior.probs, 
    evidences = c(
      Evidence_lm(y = y.in, x = x.in, model = model0, error.var = error.var), 
      Evidence_lm(y = y.in, x = x.in, model = model1, error.var = error.var)
    )
  )
  # get new data
  x.new = rep(NA, n)
  y.new = rep(NA, n)
  post.probs.mat = matrix(rep(post.probs0, n + 1), nrow = n + 1, ncol = 2, 
                          byrow = TRUE)
  post.probs.cur = post.probs0
  y.cur = y.in
  x.cur = x.in
  if(n > 0){
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
      # get y
      y.new[i] = simulateY_fromfunction(
        x = x.new[i], true.function = true.function, error.var = error.var)
      # update current information
      x.cur = c(x.cur, x.new[i])
      y.cur = c(y.cur, y.new[i])
      post.probs.cur = getHypothesesPosteriors(
        prior.probs = prior.probs, 
        evidences = c(
          Evidence_lm(y.cur, x.cur, model0, error.var),
          Evidence_lm(y.cur, x.cur, model1, error.var)
        )
      )
      post.probs.mat[i + 1, ] = post.probs.cur
    }
  }
  
  return(list(
    x.in = x.in, 
    y.in = y.in,
    x.new = x.new,
    y.new = y.new,
    post.probs = post.probs.mat
  ))
}
