# require("construct_design_matrix.R")

# for 1 and 2 dimensions

##################
### evaluate D ###
##################

### --- Posterior Variance --- ###

postvar = function(D, N, var_e, var_beta, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type)
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * (1 / var_beta)))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_beta)))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    if(diagPrior == TRUE) return(var_e * solve(crossprod(X) + var_e * diag(1/diag(var_beta))))
    else return(var_e * solve(crossprod(X) + var_e * solve(var_beta)))
  }
}

