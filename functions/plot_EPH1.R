EPH_seq_1design = function(D, N, true_beta, mean_beta0, mean_beta1, var01_seq, var_e,
                           numSims, true_model_type, H01_model_types, seed = 123){
  EPH = matrix(NA, 3, length(var01_seq))
  for(i in 1:length(var01_seq)){
    var01.temp = var01_seq[i]
    var_beta0.temp = diag(rep(var01.temp,length(mean_beta0)))
    var_beta1.temp = diag(rep(var01.temp,length(mean_beta1)))
    EPH[ , i] = calcExpPostProbH(D, N, true_beta, mean_beta0, mean_beta1, var_beta0.temp, var_beta1.temp, var_e,
                                 numSims, true_model_type, H01_model_types, seed)
  }
  rownames(EPH) = c("E[P(H0|Y,D)]", "E[P(H1|Y,D)]", "E[BF01|Y,D]")
  colnames(EPH) = paste("sigmasq01=", var01_seq, sep = "")
  return(EPH)
}

EPH_seq = function(designs, N, true_beta, mean_beta0, mean_beta1, var01_seq, var_e,
                   numSims, true_model_type, H01_model_types, seed = 123){
  # "designs" is a list of designs
  numDesigns = length(designs)
  len_var01_seq = length(var01_seq)
  # EPHs = list()
  EPH0 = matrix(NA, numDesigns, len_var01_seq)
  EPH1 = matrix(NA, numDesigns, len_var01_seq)
  EBF01 = matrix(NA, numDesigns, len_var01_seq)
  for(k in 1:numDesigns){
    D = designs[[k]]
    result = EPH_seq_1design(D, N, true_beta, mean_beta0, mean_beta1, var01_seq, var_e, 
                             numSims, true_model_type, H01_model_types, seed)
    # EPHs[[k]] = result
    EPH0[k, ] = result[1, ]
    EPH1[k, ] = result[2, ]
    EBF01[k, ] = result[3, ]
  }
  colnames(EPH0) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(EPH1) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(EBF01) = paste("sigmasq01=", var01_seq, sep = "")
  rownames(EPH0) = paste("des", 1:numDesigns, sep = "")
  rownames(EPH1) = paste("des", 1:numDesigns, sep = "")
  rownames(EBF01) = paste("des", 1:numDesigns, sep = "")
  return(list("EPH0" = EPH0, "EPH1" = EPH1, "EBF01" = EBF01))
}

plot_EPH1 = function(designs, N, betaT, mean_beta0, mean_beta1, var01_seq, var_e,
                     numSims, true_model_type, H01_model_types, seed = 123, cols = NULL){
  EPH_des = EPH_seq(designs, N, betaT, mean_beta0, mean_beta1, var01_seq, var_e,
                    numSims, true_model_type, H01_model_types, seed = 123)
  plot_ylim = range(EPH_des$EPH1[!is.nan(EPH_des$EPH1)])
  if(is.null(cols)){
    plot(x = var01_seq, y = EPH_des$EPH1[1, ], type = "l", ylim = plot_ylim, 
         ylab = "E[P(H1|Y,D)]",  xlab = "sigmasq01")
    for(i in 2:length(designs)){
      lines(x = var01_seq, y = EPH_des$EPH1[i, ], col = i)
    }
  } else{
    if(length(cols) != length(designs)) stop("length of cols argument does not match length of designs list")
    plot(x = var01_seq, y = EPH_des$EPH1[1, ], type = "l", ylim = plot_ylim, 
         ylab = "E[P(H1|Y,D)]",  xlab = "sigmasq01", col = cols[1])
    for(i in 2:length(designs)){
      lines(x = var01_seq, y = EPH_des$EPH1[i, ], col = cols[i])
    }
  }
}
