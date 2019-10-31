
constructDesignX = function(D, N, type){
  if(is.null(dim(D))){
    if(N != length(D)) stop("N isn't the same length as D")
    X = NULL
    if(type == 1) X = D
    if(type == 2) X = cbind(rep(1, N), D)
    if(type == 3) X = cbind(rep(1, N), D, D^2)
    if(type == 4) X = cbind(rep(1, N), D, D^2, D^3)
  } else{
    if(N != dim(D)[1]) stop("N is not equal to the number of rows in D")
    X = NULL
    if(type == 4){
      X = cbind(rep(1, N), D)
    }
    if(type == 5){
      X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
    }
  }
  return(X)
}
