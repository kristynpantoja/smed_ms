getGPPredictive = function(
  x, x.input, y.input, type, l, p, signal.var = 1, measurement.var = NULL){
  k_star = t(getCov(
    X1 = x, X2 = x.input, type = type, l = l, p = p, signal.var = signal.var))
  if(is.null(measurement.var)){
    K_obs = getCov(
      X1 = x.input, X2 = x.input, type = type, l = l, p = p,
      signal.var = signal.var)
  } else{
    K_obs = getCov(
      X1 = x.input, X2 = x.input, type = type, l = l, p = p, 
      signal.var = signal.var) + measurement.var * diag(length(y.input))
  }
  pred_mean = t(k_star) %*% solve(K_obs, y.input)
  pred_cov = getCov(
    X1 = x, X2 = x, type = type, l = l, p = p, signal.var = signal.var) - 
    (t(k_star) %*% solve(K_obs, k_star))
  return(list(
    "x" = x, 
    "pred_mean" = as.vector(pred_mean), 
    "pred_var" = pred_cov))
}
