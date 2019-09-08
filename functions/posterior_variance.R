##################
### evaluate D ###
##################

### --- Posterior Variance --- ###

postvar = function(D, N, var_e, var_mean, type, diagPrior = TRUE){
  X = NULL
  if(type == 1) X = D
  if(type == 2) X = cbind(rep(1, N), D)
  if(type == 3) X = cbind(rep(1, N), D, D^2)
  if(type == 4){
    X = cbind(rep(1, N), D)
  }
  if(type == 5){
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
  }
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * (1 / var_mean)))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * diag(1/diag(var_mean))))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_mean)))
  }
}

