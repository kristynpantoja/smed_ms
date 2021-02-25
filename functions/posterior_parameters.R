# works for multidimensional case too

# posterior mean of beta
postmean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       type = NULL, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% 
    (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}

# posterior variance of beta
postvar = function(D, N, var_e, var_beta, type = NULL, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    if(diagPrior == TRUE){
      return(var_e * solve(crossprod(X) + var_e * (1 / var_beta)))
    }
    else return(var_e * solve(crossprod(X) + var_e * solve(var_beta)))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    if(diagPrior == TRUE){
      return(var_e * solve(crossprod(X) + var_e * diag(1/diag(var_beta))))
    } else{
      return(var_e * solve(crossprod(X) + var_e * solve(var_beta)))
    }
  }
}
