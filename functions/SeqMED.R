# require("charge_function_q.R")
# require("construct_design_matrix.R")
# require("wasserstein_distance.R")
# require("posterior_parameters.R")
# require("add_MMED.R")
# require("simulate_y.R")

################################################################################
# simulate SeqMED, using preliminary data (or, if none, generate some)
################################################################################


# simulate_seqMED
SeqMED = function(
  D1 = NULL, y1 = NULL, Nprelim, true_beta, true_type, mean_beta0, mean_beta1, 
  var_beta0, var_beta1, var_e, f0 = NULL, f1 = NULL, type = NULL, 
  numCandidates = 10^5, k = 4, xmin = 0, xmax = 1, p = 1, 
  numSeq = 5, seqN = 10, alpha_seq = NULL, 
  buffer_seq = 0, wasserstein0 = 1, genCandidates = 1, candidates = NULL, 
  seed = NULL
  ){
  if(!is.null(seed)) set.seed(seed)
  # some checks
  if(is.null(D1)){
    if(is.null(alpha_seq)){ # generate space-filling design for first step
      D1 = MMED(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, seqN[1], 
                numCandidates, k, xmin, xmax, p, alpha = 0, buffer_seq[1], 
                genCandidates = 1, initialpt = 1)
    } else{
      D1 = MMED(mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, f0, f1, type, seqN[1], 
                numCandidates, k, xmin, xmax, p, alpha_seq[1], buffer_seq[1], 
                genCandidates = 1, initialpt = 1)
    }
  }
  if(numSeq > 1 & length(seqN) == 1) seqN = rep(seqN, numSeq)
  if(numSeq > 1 & length(buffer_seq) == 1) buffer_seq = rep(buffer_seq, numSeq)
  if(numSeq > 1 & is.null(alpha_seq)) alpha_seq = 1 # alpha_seq = c(0, ((2:numSeq) / numSeq) * (2 * p))
  if(numSeq > 1 & length(alpha_seq) == 1) alpha_seq = rep(alpha_seq, numSeq)
  
  # get y1
  if(is.null(y1)) y1 = as.vector(simulateY(D1, seqN[1], true_beta, var_e, 1, true_type))
  Nttl = sum(seqN)
  D = D1
  y = y1
  
  current_postvar0 = diag(postvar(D, length(D), var_e, var_beta0, type[1]))
  current_postmean0 = postmean(y, D, length(D), mean_beta0, var_beta0, var_e, type[1])
  current_postvar1 = diag(postvar(D, length(D), var_e, var_beta1, type[2]))
  current_postmean1 = postmean(y, D, length(D), mean_beta1, var_beta1, var_e, type[2])
  
  postvar0 = matrix(current_postvar0, length(mean_beta0), numSeq)
  postmean0 = matrix(current_postmean0, length(mean_beta0), numSeq)
  postvar1 = matrix(current_postvar1, length(mean_beta1), numSeq)
  postmean1 = matrix(current_postmean1, length(mean_beta1), numSeq)
  
  if(numSeq == 1){
    return(list("D" = D, "y" = y, 
                "postvar0" = current_postvar0, "postmean0" = current_postmean0, 
                "postvar1" = current_postvar1, "postmean1" = current_postmean1))
  }
  
  # -- Generate Candidate Points -- #
  if(is.null(candidates)){
    if(genCandidates == 1) candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
    if(genCandidates == 2) candidates = sort(runif(numCandidates, min = xmin, max = xmax))
  }
  
  print(paste("finished ", 1, " out of ", numSeq, " steps", sep = ""))
  for(t in 2:numSeq){
    
    Dt = SeqMED_batch(D, y, mean_beta0, mean_beta1, 
                      var_beta0, var_beta1, var_e,  f0, f1, type, seqN[t], 
                      numCandidates, k, xmin, xmax, p, alpha_seq[t], buffer_seq[t],
                      wasserstein0, genCandidates, candidates)
    
    yt = as.vector(simulateY(Dt$addD, seqN[t], true_beta, var_e, 1, true_type))
    
    # update D and y with new data
    D = c(D, Dt$addD)
    y = c(y, yt)
    
    # also save posterior means and variances
    postvar0[ , t] = diag(postvar(D, length(D), var_e, var_beta0, type[1]))
    postmean0[ , t] = postmean(y, D, length(D), mean_beta0, var_beta0, var_e, type[1])
    postvar1[ , t] = diag(postvar(D, length(D), var_e, var_beta1, type[2]))
    postmean1[ , t] = postmean(y, D, length(D), mean_beta1, var_beta1, var_e, type[2])
    
    print(paste("finished ", t, " out of ", numSeq, " steps", sep = ""))
  }
  return(list("D" = D, "y" = y, "postvar0" = postvar0, "postmean0" = postmean0, 
              "postvar1" = postvar1, "postmean1" = postmean1))
}