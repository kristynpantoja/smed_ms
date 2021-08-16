################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 1 or 2
#   where H1 is true

dimT = 2 # 2, 3

if_plot = FALSE

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_gpvs"
data_dir = paste0(sims_dir, "/simulated_data/outputs")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

# for generating initial data
# source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))

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
numSims = 100
numSeq = 100
seqN = 1
N = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)

# hypothesis settings
type01 = c(2, 3)
sigmasq = 0.1
mu_full = c(0.5, 0.5, 0.5) #
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# generate functions
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
