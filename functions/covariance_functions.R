### kernel functions ---

phi_sqexp = function(Xi, Xj, l){
  exp(-0.5 * (sum((Xi - Xj)^2)) / l ^ 2) # Xi, Xj could be vectors
}
# phi_sqexp = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)

# exponential, a.k.a. matern, nu = 1/2
phi_exp = function(Xi, Xj, l){
  exp(- sum(abs(Xi - Xj)) / l)
}
# matern, nu = 3/2
phi_matern2 = function(Xi, Xj, l){
  r = sum(abs(Xi - Xj))
  (1 + (sqrt(3) * r / l)) * exp(- sqrt(3) * r / l)
}
# general matern, using RandomFieldsUtils package
phi_matern = function(Xi, Xj, l, nu = 3 / 2){
  RandomFieldsUtils::matern(sum(abs(Xi - Xj)), nu, scaling = "matern")
}

# a periodic kernel that allows to adjust the period, p. 
# https://www.cs.toronto.edu/~duvenaud/cookbook/
phi_periodic = function(Xi, Xj, l, p){
  r = sum(abs(Xi - Xj))
  sinsq_arg = pi * r / p
  exp(-2 * (sin(sinsq_arg))^2 / l^2)
}

### covariance matrix ---

getCov = function(X1, X2, type, l = 1, p = 0.26, signal.var = 1){
  if(is.null(p)) p = 0.26
  phi = NULL
  if(type == 1 | type == "squaredexponential" | type == "sqexp" | 
     type == "gaussian" | type == "s" | type == "g") {
    phi = phi_sqexp
  } else if(type == 2) {
    phi = phi_matern
  } else if(type == 3 | type == "exp" | type == "exponential" | type == "e") {
    phi = phi_exp
  } else if(type == 4 | type == "matern" | type == "mat" | type == "m") {
    phi = phi_matern2
  } else if(type == 5 | type == "periodic" | type == "p") {
    phi = function(X1, X2, l){
      phi_periodic(X1, X2, l, p = p)
    }
  } else{
    stop("invalid type specification of covariance function for GP")
  }
  if(is.vector(X1) & is.vector(X2)) {
    X1 = as.matrix(X1)
    X2 = as.matrix(X2)
  }
  Sigma = matrix(NA, nrow=dim(X1)[1], ncol = dim(X2)[1])
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] = signal.var * 
        phi(X1[i, , drop = FALSE], X2[j, , drop = FALSE], l)
    }
  }
  return(Sigma)
}
