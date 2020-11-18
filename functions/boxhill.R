# D, criterion in Box and Hill 1967
#   for y ~ N(eta, sigma^2)
#     eta = E[y] = function that is linear (or approx linear) in parameters
#   and for m models

# get posterior distribution for beta
getBetaPosterior = function(y, X, beta_prior_mean, beta_prior_var, var_e, diagPrior = TRUE){
  # y : response
  # X : design matrix
  
  # calculate the posterior variance (a matrix)
  beta.var.inv = NA
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
  postmean = (1 / var_e) * postvar %*% 
    (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(list(mean = postmean, var = postvar))
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

# for m = 2 models, i and j:
BHD_m2 = function(
  y, # observed responses
  post.probs, # 2-vector
  error.var, # known error variance
  model.i, # design matrix of covariates for observed points X, 
  #   covariates of new point X.n, beta.mean, beta.var -- for model i
  model.j # design matrix of covariates for observed points X, 
  #   covariates of new point X.n, beta.mean, beta.var -- for model j
){
  # posterior predictive distributions
  post.i = getBetaPosterior( # beta's posterior, under model i
    y, model.i$X, model.i$beta.mean, model.i$beta.var, error.var
    )
  postpred.i = list( # posterior predictive
    mean = model.i$X.n %*% post.i$mean, 
    var = diag(post.i$var) + error.var
  )
  post.j = getBetaPosterior( # beta's posterior, under model i
    y, model.j$X, model.j$beta.mean, model.j$beta.var, error.var
    )
  postpred.j = list( # posterior predictive
    mean = model.j$X.n %*% post.j$mean, 
    var = diag(post.j$var) + error.var
  )
  # evaluate criterion D
  inner.sigma.term = (postpred.i$var - postpred.j$var)^2 / 
    ((error.var + postpred.i$var) * (error.var + postpred.j$var))
  inner.y.term = (postpred.i$mean - postpred.j$mean)^2 * # diff bt means
    ((error.var + postpred.i$var)^(-1) + (error.var + postpred.j$var)^(-1))
  bhd = 0.5 * prod(post.probs) * (inner.sigma.term + inner.y.term)
  return(bhd)
}

# BH_m2_add1 = function(
#   y, # observed responses
#   prior.probs, # 2-vector
#   error.var, # known error variance
#   model.i, # design matrix of covariates for observed points X, 
#   #   covariates of new point X.n, beta.mean, beta.var -- for model i
#   model.j # design matrix of covariates for observed points X, 
#   #   covariates of new point X.n, beta.mean, beta.var -- for model j
# ){
#   # make sure number of elements in initD matches number of elements in y
#   if(length(y) != dim(model.i$X)[1] | length(y) != dim(model.j$X)[1]){
#     stop("Error in BHgp_m2_add1() : length of y input doesn't match 
#          length of x input from one of the pairs of models!!")
#   }
#   
#   # evaluate criterion over x_seq
#   bhd_seq = BHD_m2(y, post.probs, error.var, model.i, model.j)
#   if(!all(!is.nan(bhd_seq))) warning("Warning in BHD_m2() : There were NaNs in Box & Hill criterion evaluation over the candidate set!!")
#   # which(is.nan(BHD_seq))
#   
#   # get new point
#   x.new.idx = which.max(bhd_seq)
#   x.new = x_seq[x.new.idx]
#   
#   return(list(
#     x.new = x.new,
#     x.new.idx = x.new.idx
#   ))
# }

BH_m2 = function(
  y, # observed responses
  prior.probs, # 2-vector
  error.var, # known error variance
  model.i, # design matrix of covariates for observed points X, 
  #   covariates of new point X.n, beta.mean, beta.var -- for model i
  model.j # design matrix of covariates for observed points X, 
  #   covariates of new point X.n, beta.mean, beta.var -- for model j
){
  
  # make sure number of elements in initD matches number of elements in y
  if(length(y) != dim(model.i$X)[1] | length(y) != dim(model.j$X)[1]){
    stop("Error in BHgp_m2_add1() : length of y input doesn't match 
         length of x input from one of the pairs of models!!")
  }
  
  # posterior predictive distributions
  post.i = getBetaPosterior( # beta's posterior, under model i
    y, model.i$X, model.i$beta.mean, model.i$beta.var, error.var
  )
  postpred.i = list( # posterior predictive
    mean = model.i$X.n %*% post.i$mean, 
    var = diag(post.i$var) + error.var
  )
  post.j = getBetaPosterior( # beta's posterior, under model i
    y, model.j$X, model.j$beta.mean, model.j$beta.var, error.var
  )
  postpred.j = list( # posterior predictive
    mean = model.j$X.n %*% post.j$mean, 
    var = diag(post.j$var) + error.var
  )
  # evidences: p_i(nth observation), where p_i is the density under model i
  evidences = c(
    dnorm(y, postpred.i$mean, postpred.i$var), 
    dnorm(y, postpred.j$mean, postpred.j$var)
  )
  # posterior probabilities of H0, H1 using preliminary data
  post.probs0 = getHypothesesPosteriors( # posterior probability with current data
    prior.probs = prior.probs, 
    evidences = evidences
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
    bhd_seq = BHD_m2(y.cur, post.probs.cur, error.var, model0, model1)
    if(!all(!is.nan(bhd_seq))){
      warning("Warning in BHDgp_m2() : There were NaNs in Box & Hill criterion 
              evaluation over the candidate set!!")
    }
    # get new point
    x.new.idx[i] = which.max(bhd_seq)
    x.new[i] = x_seq[x.new.idx[i]]
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
  return(BHD)
}