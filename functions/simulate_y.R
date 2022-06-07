
simulateY_fromfunction = function(
  x, true.function, error.var, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  Y = true.function(x) + rnorm(length(x), 0, sqrt(error.var))
  return(as.vector(Y))
}

simulateY_frommultivarfunction = function(
  x, true.function, error.var, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  # if x isn't in matrix form, turn it into a matrix
  if(!("matrix" %in% class(x))) x = t(matrix(x))
  Y = true.function(x) + rnorm(nrow(x), 0, sqrt(error.var))
  return(as.vector(Y))
}
