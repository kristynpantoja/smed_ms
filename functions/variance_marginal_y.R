# Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly = function(x, var_mean, var_e, type, var_margy){
  # type:
  #   1 for linear model without slope
  #   2 for linear model with slope
  #   3 for quadratic model with slope
  if(!is.null(type)){
    if(type == 1) var_e + x^2 * var_mean
    else if(type == 2){
      vars = diag(var_mean)
      vars[1] + x^2 * vars[2] + var_e
    }
    else if(type == 3){
      vars = diag(var_mean)
      vars[1] + x^2 * vars[2] + x^4 * vars[3] + var_e
    }
    else stop(paste("invalid type given : ", type))
  } else{
    var_margy(x = x, var_mean = var_mean, var_e = var_e)
  }
  
}

# Var[y | H_m], after marginalizing out \beta, for some hypothesis m
var_marginaly_2d = function(x, var_mean, var_e, type, var_margy){
  # type:
  #   1 for linear model without slope
  #   2 for linear model with slope
  #   3 for quadratic model with slope
  #   4 for two-dimentional model
  if(!is.null(type)){
    if(type == 4){
      vars = diag(var_mean)
      vars[1] + x[1]^2 * vars[2] + x[2]^2 * vars[3] + var_e
    }
    else if(type == 5){
      vars = diag(var_mean)
      vars[1] + x[1]^2 * vars[2] + x[1]^4 * vars[3] + x[2]^2 * vars[4] + x[2]^4 * vars[5] + var_e
    }
    else stop(paste("invalid type given : ", type))
  } else{
    if(!is.null(var_margy)) var_margy(x = x, var_mean = var_mean, var_e = var_e)
    else stop("cannot compute var_marginaly: no type given, and no var_margy fn given")
  }
}