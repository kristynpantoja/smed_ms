# for 1 and 2 dimensions

##################
### evaluate D ###
##################

### --- Posterior Variance --- ###

postvar = function(D, N, var_e, var_mean, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * (1 / var_mean)))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * diag(1/diag(var_mean))))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  }
}

