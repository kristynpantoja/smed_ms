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

