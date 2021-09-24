SeqMEDvs = function(
  y.in = NULL, x.in = NULL, model0, model1, error.var, 
  candidates, true.function, true.indices, dimX, 
  k = 4, xmin = -1, xmax = 1, p = 1, 
  numSeq = 5, seqN = 10, alpha_seq = 1, 
  prints = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
  if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  if(is.null(candidates)){
    if(genCandidates == 1){
      candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    } else if(genCandidates == 2){
      candidates = sort(runif(numCandidates, min = xmin, max = xmax))
    } else{
      stop("SeqMED: invalid genCandidates argument")
    }
  }
  
  x.cur = x.in
  y.cur = y.in
  x.new = matrix(NA, nrow = 0, ncol = dimX)
  y.new = c()
  
  if(numSeq == 1){
    # get posterior distributions of beta
    postbeta0 = getBetaPosterior(
      y = y.cur, X = x.cur[, model0$indices, drop = FALSE], model0$beta.mean, 
      model0$beta.var, error.var)
    postbeta1 = getBetaPosterior(
      y = y.cur, X = x.cur[, model1$indices, drop = FALSE], model1$beta.mean, 
      model1$beta.var, error.var)
    
    return(list(
      x.in = x.in, 
      y.in = y.in, 
      x.new = x.new,
      y.new = y.new,
      postvar0 = diag(postbeta0$var), postmean0 = postbeta0$mean, 
      postvar1 = diag(postbeta1$var), postmean1 = postbeta1$mean))
  }
  
  for(t in 1:numSeq){
    Dt = SeqMEDvs_batch(
      initD = x.cur, inity = y.cur, model0 = model0, model1 = model1, 
      error.var = error.var, N2 = seqN[t], 
      candidates = candidates, dimX = dimX, xmin = xmin, xmax = xmax, 
      k = k, p = p, alpha = alpha_seq[t], 
      batch.idx = t)
    
    yt = simulateY_frommultivarfunction(
      x = Dt$addD, true.function = true.function, 
      error.var = error.var)
    
    # update x.cur and y.cur with new data
    x.cur = rbind(x.cur, Dt$addD)
    y.cur = c(y.cur, yt)
    x.new = rbind(x.new, Dt$addD)
    y.new = c(y.new, yt)
    
    if(prints){
      print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
    }
  }
  
  # get updated posterior distributions of beta
  #   here, x.cur is rbind(x.in, x.new) & y.cur = c(x.in, y.new)
  #   i.e. the cumulative design
  postbeta0 = getBetaPosterior(
    y = y.cur, X = x.cur[, model0$indices, drop = FALSE], model0$beta.mean, 
    model0$beta.var, error.var)
  postbeta1 = getBetaPosterior(
    y = y.cur, X = x.cur[, model1$indices, drop = FALSE], model1$beta.mean, 
    model1$beta.var, error.var)
  
  return(list(
    x.in = x.in, 
    y.in = y.in, 
    x.new = x.new,
    y.new = y.new,
    postvar0 = diag(postbeta0$var), postmean0 = postbeta0$mean, 
    postvar1 = diag(postbeta1$var), postmean1 = postbeta1$mean))
}


















































































































































































































# 
# # for multidimensional case... like LASSO or ridge
# generate_SMMEDvs = function(mean_beta_full, beta_true = NULL, indices_true = NULL, 
#                             indices0, indices1, mean_beta0 = NULL, mean_beta1 = NULL, 
#                             var_e = 1, var_beta = NULL, var_beta0 = NULL, var_beta1 = NULL, 
#                             xmin = -1, xmax = 1, numCandidates = 10^5, k = 4, p = 1,
#                             initD = NULL, inity = NULL, numSeq = 5, N_seq = 10, 
#                             alpha_seq = NULL, buffer_seq = 0, candidates = NULL, 
#                             wasserstein0 = 1, genCandidates = 1, seed = NULL){
#   
#   # first check if some things are null
#   # some checks!!!!!!!!!!!!!!!!!!!!!!
#   if(is.null(mean_beta_full)) stop("mean_beta_full is not given")
#   if(is.null(beta_true)){
#     if(is.null(indices_true)){
#       stop("neither beta_true nor indices_true are given")
#     } else{
#       beta_true = mean_beta_full[indices_true]
#     }
#   } else{
#     if(!is.null(indices_true)){
#       if(!all.equal(beta_true, mean_beta_full[indices_true])){
#         stop("beta_true does not match indices_true")
#       }
#     }
#   }
#   p0 = length(mean_beta0)
#   p1 = length(mean_beta1)
#   if(is.null(var_beta0) & is.null(var_beta1)){
#     if(is.null(var_beta)){
#       stop("neither var_beta0, var_beta1, nor var_beta are given")
#     } else{
#       var_beta0 = diag(rep(var_beta, p0))
#       var_beta1 = diag(rep(var_beta, p1))
#     }
#   }
#   if(is.null(initD)){
#     stop("no initial data! haven't written this case for multiple dimensions yet. going to just randomly select points.")
#     # if(is.null(alpha_seq)){ # generate space-filling design for first step
#     #   D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
#     #                          numCandidates, k, xmin, xmax, p, alpha = 0, buffer_seq[1], 
#     #                          genCandidates = 1, initialpt = 1)
#     # } else{
#     #   D1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, N_seq[1], 
#     #                          numCandidates, k, xmin, xmax, p, alpha_seq[1], buffer_seq[1], 
#     #                          genCandidates = 1, initialpt = 1)
#     # }
#     # initN = length(D1)
#   } else{
#     D1 = initD
#     initN = dim(D1)[1]
#   }
#   # sequence checks
#   if(length(N_seq) == 1) N_seq = rep(N_seq, numSeq)
#   if(length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
#   if(is.null(alpha_seq)) alpha_seq = rep(1, numSeq)
#   if(!is.null(alpha_seq) & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
#   
#   if(!is.null(seed)) set.seed(seed)
#   # get y1 (if not given)
#   if(is.null(inity)){
#     y1 = as.vector(simulateYvs(D1[ , indices_true], N_seq[1], beta_true, sigmasq, 1, seed = seed))
#   } else{
#     y1 = inity
#   }
#   y = y1
#   D = D1
#   if(numSeq == 1){
#     return(list("D" = D1,
#                 "y" = y1))
#   }
#   
#   # why do I do this? why only take the diagonal of postvar?
#   current_postvar0 = diag(postvar(D[ , indices0], initN, var_e, var_beta0, type = NULL))
#   current_postmean0 = postmean(y, D[ , indices0], initN, mean_beta0, var_beta0, var_e, type = NULL)
#   current_postvar1 = diag(postvar(D[ , indices1], initN, var_e, var_beta1, type = NULL))
#   current_postmean1 = postmean(y, D[ , indices1], initN, mean_beta1, var_beta1, var_e, type = NULL)
#   
#   # to save it more compactly; only need pointwise variance, no need to save whole covariance matrix
#   # first one is for initial data
#   # the other numSeq - 1 ones will be replaced during the numSeq loop
#   postvar0 = matrix(current_postvar0, length(mean_beta0), numSeq)
#   postmean0 = matrix(current_postmean0, length(mean_beta0), numSeq)
#   postvar1 = matrix(current_postvar1, length(mean_beta1), numSeq)
#   postmean1 = matrix(current_postmean1, length(mean_beta1), numSeq)
#   
#   # get candidates
#   dimX = length(mean_beta_full)
#   if(is.null(candidates)){
#     candidates = matrix(runif(n = dimX * numCandidates, min = xmin, max = xmax), 
#                         nrow = numCandidates, ncol = dimX)
#   }
#   
#   print(paste("finished ", 0, " out of ", numSeq, " steps", sep = ""))
#   for(t in 1:numSeq){
#     seed = seed + t
#     
#     Dt = SeqMEDvs_batch(
#       initD = D, inity = y, mean_beta_full = mean_beta_full, 
#       beta_true = beta_true, indices_true = indices_true, 
#       indices0 = indices0, indices1 = indices1, 
#       mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
#       var_e = var_e, var_beta = var_beta, 
#       var_beta0 = var_beta0, var_beta1 = var_beta1,
#       N2 = N_seq[t], xmin = xmin, xmax = xmax, 
#       numCandidates = numCandidates, k = k, p = p, 
#       alpha = alpha_seq[t], buffer = buffer_seq[t], candidates,
#       wasserstein0 = wasserstein0, genCandidates = genCandidates, 
#       dimX = dimX)
#     
#     yt = as.vector(simulateYvs(Dt$addD[ , indices_true, drop = FALSE], N_seq[t], beta_true, sigmasq, 1, seed = seed))
#     
#     # update D and y with new data
#     D = rbind(D, Dt$addD)
#     y = c(y, yt)
#     
#     # also save posterior means and variances
#     currentN = dim(D)[1]
#     postvar0[ , t] = diag(postvar(D[ , indices0], currentN, var_e, var_beta0, type = NULL))
#     postmean0[ , t] = postmean(y, D[ , indices0], currentN, mean_beta0, var_beta0, var_e, type = NULL)
#     postvar1[ , t] = diag(postvar(D[ , indices1], currentN, var_e, var_beta1, type = NULL))
#     postmean1[ , t] = postmean(y, D[ , indices1], currentN, mean_beta1, var_beta1, var_e, type = NULL)
#     
#     print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
#   }
#   return(list("initD" = initD, "addD" = D[(initN + 1):dim(D)[1], ], "D" = D, "y" = y, 
#               "postvar0" = postvar0, "postmean0" = postmean0, 
#               "postvar1" = postvar1, "postmean1" = postmean1))
# }
# 
# 
