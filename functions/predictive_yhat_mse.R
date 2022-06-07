getMSEYhat = function(
  pt, D, N, beta.true, true_type, beta.mean, beta.var, error.var, 
  prior_type, indices = NULL, diagPrior = TRUE
){
  # note true_type must match prior_type
  if(!is.null(true_type)) if(true_type != prior_type) warning("true_type is not the same as prior_type: will lead to non-conforming matrices in multiplication")
  if(!is.null(true_type)){
    # variance term
    x.temp = matrix(constructDesignX(pt, 1, prior_type))
    X.temp = constructDesignX(D, N, prior_type) # design matrix, depends on the model (type)
    x.true = matrix(constructDesignX(pt, 1, true_type))
    Sigma_B = postvar(D, N, error.var, beta.var, prior_type) # posterior variance
    XtX = crossprod(X.temp)
  } else{
    # variance term
    x.temp = matrix(constructDesignX(pt, 1, prior_type))[indices]
    X.temp = constructDesignX(D, N, prior_type)[ , indices]
    x.true = matrix(constructDesignX(pt, 1, true_type))[indices]
    Sigma_B = postvar(D[ , indices], N, error.var, beta.var, prior_type)[ , indices] # posterior variance
    XtX = crossprod(X.temp)
  }
  varBn = (1/error.var) * Sigma_B %*% XtX %*% Sigma_B
  variance_term = t(x.temp) %*% varBn %*% x.temp
  
  # bias-squared term
  # calculate expectation of posterior mean
  expect_postmean = (1/error.var) * Sigma_B %*% XtX %*% matrix(beta.true) + 
    Sigma_B %*% solve(beta.var) %*% matrix(beta.mean)
  biassq_term = (t(x.temp) %*% expect_postmean - t(x.true) %*% beta.true)^2
  
  # MSE of y-hat
  MSEyhat = variance_term + biassq_term
  return(list("var_term" = variance_term, "biassq_term" = biassq_term, "MSEyhat" = as.vector(MSEyhat)))
}

getMSEYhat_seq = function(
  x_seq, D, N, beta.true, true_type, beta.mean, beta.var, error.var, 
  prior_type, indices = NULL, diagPrior = TRUE
){
  mseyhat_seq = rep(NA, length(x_seq))
  mseyhat_var_seq = rep(NA, length(x_seq))
  mseyhat_biassq_seq = rep(NA, length(x_seq))
  for(i in 1:length(x_seq)){
    mseyhat = getMSEYhat(
      x_seq[i], D, N, beta.true, true_type, beta.mean, beta.var, 
      error.var, prior_type, indices, diagPrior)
    mseyhat_seq[i] = mseyhat$MSEyhat
    mseyhat_var_seq[i] = mseyhat$var_term
    mseyhat_biassq_seq[i] = mseyhat$biassq_term
  }
  return(list(
    "MSEyhat" = mseyhat_seq, 
    "var_term" = mseyhat_var_seq, 
    "biassq_term" = mseyhat_biassq_seq
  ))
}
