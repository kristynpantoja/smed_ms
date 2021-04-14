getGPPredictive = function(
  x, x.input, y.input, type, l, signal.var = 1, error.var = NULL){
  k_star = t(getCov(x, x.input, type, l, signal.var))
  if(is.null(error.var)){
    K_obs = getCov(x.input, x.input, type, l, signal.var)
  } else{
    K_obs = getCov(x.input, x.input, type, l, signal.var) + 
      error.var * diag(length(y.input))
  }
  pred_mean = t(k_star) %*% solve(K_obs, y.input)
  pred_cov = getCov(x, x, type, l) - (t(k_star) %*% solve(K_obs, k_star))
  return(list(
    "x" = x, 
    "pred_mean" = as.vector(pred_mean), 
    "pred_var" = pred_cov))
}