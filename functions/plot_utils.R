####################################################
# plot EPPH ########################################
####################################################

# require("construct_design_matrix.R")
# require("simulate_y.R")

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


####################################################
# plot MSE #########################################
####################################################


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


####################################################
# plot posterior variance ##########################
####################################################

postvar_seq_1design = function(D, N, mean_beta, var01_seq, var_e, model_type){ # doesn't depend on true beta
  p = length(mean_beta) # note: model type isn't necessarily true model type
  postvars = matrix(NA, p, length(var01_seq))
  for(i in 1:length(var01_seq)){
    var01.temp = var01_seq[i]
    var_beta.temp = diag(rep(var01.temp,length(mean_beta)))
    postvars[ , i] = tryCatch({diag(postvar(D, N, var_e, var_beta.temp, model_type, diagPrior = TRUE))}, 
                              error = function(e) {rep(NA, p)})
  }
  rownames(postvars) = paste("V[B", 0:(p - 1), "|Y,D]", sep = "")
  colnames(postvars) = paste("sigmasq01=", var01_seq, sep = "")
  return(postvars)
}

postvar_seq = function(designs, N, mean_beta, var01_seq, var_e, model_type){
  # "designs" is a list of designs
  numDesigns = length(designs)
  len_var01_seq = length(var01_seq)
  # EPHs = list()
  postvarB0 = matrix(NA, numDesigns, len_var01_seq)
  postvarB1 = matrix(NA, numDesigns, len_var01_seq)
  postvarB2 = matrix(NA, numDesigns, len_var01_seq)
  postvarB3 = matrix(NA, numDesigns, len_var01_seq)
  postvarB4 = matrix(NA, numDesigns, len_var01_seq)
  for(k in 1:numDesigns){
    D = designs[[k]]
    result = postvar_seq_1design(D, N, mean_beta, var01_seq, var_e, model_type)
    # EPHs[[k]] = result
    if(model_type >= 2){
      postvarB0[k, ] = result[1, ]
      postvarB1[k, ] = result[2, ]
    }
    if(model_type >= 3){
      postvarB2[k, ] = result[3, ]
    }
    if(model_type == 5){
      postvarB3[k, ] = result[4, ]
      postvarB4[k, ] = result[5, ]
    }
  }
  colnames(postvarB0) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(postvarB1) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(postvarB2) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(postvarB3) = paste("sigmasq01=", var01_seq, sep = "")
  colnames(postvarB4) = paste("sigmasq01=", var01_seq, sep = "")
  
  rownames(postvarB0) = paste("des", 1:numDesigns, sep = "")
  rownames(postvarB1) = paste("des", 1:numDesigns, sep = "")
  rownames(postvarB2) = paste("des", 1:numDesigns, sep = "")
  rownames(postvarB3) = paste("des", 1:numDesigns, sep = "")
  rownames(postvarB4) = paste("des", 1:numDesigns, sep = "")
  return(list("postvarB0" = postvarB0, "postvarB1" = postvarB1, "postvarB2" = postvarB2,
              "postvarB3" = postvarB3, "postvarB4" = postvarB4))
}

plot_postvar = function(designs, N, mean_beta, var01_seq, var_e, model_type, which.plots = NULL, cols = NULL){
  postvar_seq = postvar_seq(designs, N, mean_beta, var01_seq, var_e, model_type)
  if(is.null(which.plots)) which.plots = 1:length(mean_beta)
  for(k in 1:length(which.plots)){
    plot_ylim = range(na.omit(postvar_seq[[which.plots[k]]][!is.nan(postvar_seq[[which.plots[k]]])]))
    if(is.null(cols)){
      plot(x = var01_seq, y = postvar_seq[[which.plots[k]]][1, ], type = "l", ylim = plot_ylim, 
           ylab = paste("V[B", which.plots[k] - 1, "|Y,D]", sep = ""),  xlab = "sigmasq01")
      for(i in 2:length(designs)){
        lines(x = var01_seq, y = postvar_seq[[which.plots[k]]][i, ], col = i)
      }
    } else{
      if(length(cols) != length(designs)) stop("length of cols argument does not match length of designs list")
      plot(x = var01_seq, y = postvar_seq[[which.plots[k]]][1, ], type = "l", ylim = plot_ylim, 
           ylab = paste("V[B", which.plots[k] - 1, "|Y,D]", sep = ""),  xlab = "sigmasq01", col = cols[1])
      for(i in 2:length(designs)){
        lines(x = var01_seq, y = postvar_seq[[which.plots[k]]][i, ], col = cols[i])
      }
    }
  }
}



