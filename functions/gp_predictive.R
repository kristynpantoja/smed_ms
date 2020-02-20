getGPPredictive = function(x_seq, x_train, y_train, type_arg, l_arg, nugget = NULL){
  k_star = t(getCov(x_seq, x_train, type_arg, l_arg))
  if(is.null(nugget)) K_obs = getCov(x_train, x_train, type_arg, l_arg)
  else K_obs = getCov(x_train, x_train, type_arg, l_arg) + diag(rep(nugget, length(x_train)))
  pred_mean = t(k_star) %*% solve(K_obs, y_train)
  pred_cov = getCov(x_seq, x_seq, type_arg, l_arg) - (t(k_star) %*% solve(K_obs, k_star))
  return(list("pred_mean" = as.vector(pred_mean), "pred_var" = pred_cov))
}