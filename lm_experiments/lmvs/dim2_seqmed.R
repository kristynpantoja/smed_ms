################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 1 or 2
#   where H1 is true

dimT = 2 # 2, 3

if_plot = FALSE

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lmvs/outputs"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
# source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
# source(paste(functions_dir, "/boxhill_vs.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

library(mvtnorm)

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

library(ggplot2)
library(reshape2)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 25
numSeq = 9
seqN = 1
N = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 1000

# other settings
sigmasq = 0.3
p = 3
k = 4 * p
initN = 5
sigmasq = 0.1

# hypothesis settings
type01 = c(2, 3)
mu_full = c(0.5, 0.5, 0.5) #
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
sigmasq01 = 0.5
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# scenario: 2 out of 3 dimensions
################################################################################
betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]

# seqmed settings

# generate seqmeds
registerDoRNG(rng.seed)
seqmed_list = foreach(i = 1:numSims) %dorng% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  SeqMED(
    D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
    beta.mean0 = mu0, beta.mean1 = mu1, beta.var0 = V0, beta.var1 = V1, 
    error.var = sigmasq, f0 = f0, f1 = f1, type = type01, xmin = xmin, xmax = xmax, 
    candidates = candidates, numSeq = numSeq, seqN = seqN
  )
  
  
  
  initD = matrix(runif(n = length(mu_full) * initN, min = xmin, max = xmax), 
                 nrow = initN, ncol = length(mu_full))
  inity = as.vector(simulateYvs(
    initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
  smmed_H0 = generate_SMMEDvs(
    mu_full, betaT, indicesT, indices0, indices1, mu0, mu1,
    sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
    initD, inity, numSeq, seqN)
  
}
saveRDS(seqmed_list, paste(output_home, "/seqmed/dim2_seqmed_simulations", 
                           "_numSeq", numSeq, 
                           "_seqN", seqN,
                           "_numSims", numSims,
                           ".rds", sep = ""))


initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
smmed_H0 = generate_SMMEDvs(mu_full, betaT, indicesT, indices0, indices1, mu0, mu1,
                            sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
                            initD, inity, numSeq, N_seq, seed = 12)


registerDoRNG(rng.seed)

if(dimT == 2){
  null_cov = getCov(x_grid, x_grid, typeT, lT, pT, sigmasq_signal) ######################## check
  null_mean = rep(0, numx^2)
} else if(dimT == 3){
  null_cov = getCov(x_seq, x_seq, typeT, lT, pT, sigmasq_signal)
  null_mean = rep(0, numx)
}

# the function values
filename_append = ""
if(is.null(sigmasq_measuremt)){
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov))
  filename_append = ""
} else{
  if(dimT == 1){
    y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov + 
                            sigmasq_measuremt * diag(numx))) 
  } else if(dimT == 2){
    y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov + 
                            sigmasq_measuremt * diag(numx^2))) 
  }
  filename_append = "_noise"
  
  # plot
  if(if_plot){
    sim_index = 3
    fields::quilt.plot(x_grid, y_seq_mat[ , sim_index])
    plot(x_grid[ , 1], y_seq_mat[ , sim_index])
  }
}

if(dimT == 1){ # expand to 2 dims
  y_seq_mat_1d = y_seq_mat
  y_seq_mat = matrix(NA, nrow = numx^2, ncol = numSims)
  for(i in 1:numSims){
    expanded = expand.grid(y_seq_mat_1d[ , i], y_seq_mat_1d[ , i])
    y_seq_mat[ , i] = expanded[ , 1]
  }
  
  # plot
  if(if_plot){
    sim_index = 3
    fields::quilt.plot(x_grid, y_seq_mat[ , sim_index])
    plot(x_grid[ , 1], y_seq_mat[ , sim_index])
    order_truef = order(x_grid[,1])
    lines(x_grid[order_truef,1], y_seq_mat[order_truef,sim_index])
  }
}

# save the results!
simulated_data_file = paste0(
  data_dir,
  "/", typeT,
  "_l", lT,
  "_dim", dimT, 
  filename_append, 
  "_seed", rng.seed,
  ".rds")
saveRDS(
  list(
    numx = numx,
    x_seq = x_seq, 
    x_grid = x_grid,
    null_mean = null_mean, 
    null_cov = null_cov, 
    numSims = numSims, 
    function_values_mat = y_seq_mat
  ), 
  file = simulated_data_file)
