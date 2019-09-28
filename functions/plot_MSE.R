MSE_seq_1design = function(D, N, true_beta, mean_beta, var01_seq, var_e, true_model_type){
  p = length(true_beta)
  MSE = matrix(NA, p, length(var01_seq))
  for(i in 1:length(var01_seq)){
    var01.temp = var01_seq[i]
    var_beta.temp = diag(rep(var01.temp,length(mean_beta)))
    MSE[ , i] = tryCatch({getClosedMSE(D, N, true_beta, mean_beta, var_beta.temp, var_e,
                                       true_model_type)$MSE_postmean}, error = function(e) {rep(NA, p)})
  }
  rownames(MSE) = paste("MSE(B", 0:(p - 1), ")", sep = "")
  colnames(MSE) = paste("sigmasq01=", var01_seq, sep = "")
  return(MSE)
}

MSE_seq = function(designs, N, true_beta, mean_beta, var01_seq, var_e, true_model_type){
  # "designs" is a list of designs
  numDesigns = length(designs)
  len_var01_seq = length(var01_seq)
  # EPHs = list()
  MSEB0 = matrix(NA, numDesigns, len_var01_seq)
  MSEB1 = matrix(NA, numDesigns, len_var01_seq)
  MSEB2 = matrix(NA, numDesigns, len_var01_seq)
  MSEB3 = matrix(NA, numDesigns, len_var01_seq)
  MSEB4 = matrix(NA, numDesigns, len_var01_seq)
  for(k in 1:numDesigns){
    D = designs[[k]]
    result = MSE_seq_1design(D, N, true_beta, mean_beta, var01_seq, var_e, true_model_type)
    # EPHs[[k]] = result
    if(true_model_type >= 2){
      MSEB0[k, ] = result[1, ]
      MSEB1[k, ] = result[2, ]
    }
    if(true_model_type >= 3){
      MSEB2[k, ] = result[3, ]
    }
    if(true_model_type == 5){
      MSEB3[k, ] = result[4, ]
      MSEB4[k, ] = result[5, ]
    }
  }
  colnames(MSEB0) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(MSEB1) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(MSEB2) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(MSEB3) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(MSEB4) = paste("sigmasq01=", var01_seq, sep = "")
  
  rownames(MSEB0) = paste("des", 1:numDesigns, sep = "")
  rownames(MSEB1) = paste("des", 1:numDesigns, sep = "")
  rownames(MSEB2) = paste("des", 1:numDesigns, sep = "")
  rownames(MSEB3) = paste("des", 1:numDesigns, sep = "")
  rownames(MSEB4) = paste("des", 1:numDesigns, sep = "")
  return(list("MSEB0" = MSEB0, "MSEB1" = MSEB1, "MSEB2" = MSEB2,
              "MSEB3" = MSEB3, "MSEB4" = MSEB4))
}

plot_MSE = function(designs, N, true_beta, mean_beta, var01_seq, var_e, true_model_type, which.plots = NULL, cols = NULL){
  MSE_des = MSE_seq(designs, N, true_beta, mean_beta, var01_seq, var_e, true_model_type)
  if(is.null(which.plots)) which.plots = 1:length(true_beta)
  for(k in 1:length(which.plots)){
    plot_ylim = range(na.omit(MSE_des[[which.plots[k]]][!is.nan(MSE_des[[which.plots[k]]])]))
    if(is.null(cols)){
      plot(x = var01_seq, y = MSE_des[[which.plots[k]]][1, ], type = "l", ylim = plot_ylim, 
           ylab = paste("MSE(B", which.plots[k] - 1, ")", sep = ""),  xlab = "sigmasq01")
      for(i in 2:length(designs)){
        lines(x = var01_seq, y = MSE_des[[which.plots[k]]][i, ], col = i)
      }
    } else{
      if(length(cols) != length(designs)) stop("length of cols argument does not match length of designs list")
      plot(x = var01_seq, y = MSE_des[[which.plots[k]]][1, ], type = "l", ylim = plot_ylim, 
           ylab = paste("MSE(B", which.plots[k] - 1, ")", sep = ""),  xlab = "sigmasq01", col = cols[1])
      for(i in 2:length(designs)){
        lines(x = var01_seq, y = MSE_des[[which.plots[k]]][i, ], col = cols[i])
      }
    }
  }
}

