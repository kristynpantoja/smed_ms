# Covariance Functions

###############################
### 1d Covariance functions ###
###############################
# radial aka squared exponential aka gaussian

phi_sqexp = function(Xi, Xj, l = 1) exp(-0.5 * (sum((Xi - Xj)^2)) / l ^ 2) # Xi, Xj could be vectors
# phi_sqexp = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)

# exponential, a.k.a. matern, nu = 1/2
phi_exp = function(Xi, Xj, l) exp(- abs(Xi - Xj) / l)
# matern, nu = 3/2
phi_matern2 = function(Xi, Xj, l){
  r = abs(Xi - Xj)
  (1 + (sqrt(3) * r / l)) * exp(- sqrt(3) * r / l)
}
# general matern, l = 1
phi_matern = function(Xi, Xj, l, nu = 3/2) RandomFieldsUtils::matern(abs(Xi - Xj), nu, scaling = "matern")
# periodic, but a generic version
# phi_periodic = function(Xi, Xj, l){
#   r = Xi - Xj
#   exp(-2 * (sin(r / 2) / l)^2)
# }
# a periodic kernel that allows to adjust the period, p. 
phi_periodic = function(Xi, Xj, l, p = pi / 24){
  r = abs(Xi - Xj)
  sinsq_arg = pi * r / p
  exp(-2 * (sin(sinsq_arg / 2) / l)^2)
}

# compute covariance matrix
# getCov = function(X, Y, type = 1, l = 1){ # change l = 1 to a list theta = NULL
#   phi = NULL
#   if(type == 1) phi = phi_sqexp
#   else if(type == 2) phi = phi_matern
#   else if(type == 3) phi = phi_exp
#   else if(type == 4) phi = phi_matern2
#   else if(type == 5) phi = phi_periodic
#   else stop("invalid type specification of covariance function for GP")
#   # outer(X, Y, FUN = Vectorize(phi), l)
#   outer(X, Y, FUN = phi, l)
# }

getCov = function(X1, X2, type = 1, l = 1){
  phi = NULL
  if(type == 1) phi = phi_sqexp
  else if(type == 2) phi = phi_matern
  else if(type == 3) phi = phi_exp
  else if(type == 4) phi = phi_matern2
  else if(type == 5) phi = phi_periodic
  else stop("invalid type specification of covariance function for GP")
  if(is.vector(X1) & is.vector(X2)) {
    X1 = as.matrix(X1)
    X2 = as.matrix(X2)
  }
  Sigma = matrix(NA, nrow=dim(X1)[1], ncol = dim(X2)[1])
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] = phi(X1[i, , drop = FALSE], X2[j, , drop = FALSE], l)
    }
  }
  return(Sigma)
}
