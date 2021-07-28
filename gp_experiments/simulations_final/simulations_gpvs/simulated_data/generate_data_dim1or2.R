################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 1 or 2
#   where H1 is true

dimT = 2
lT = 0.05
typeT = "squaredexponential"
pT = NULL

# if 1, assume 1st dimension is true.
if_plot = FALSE

################################################################################
# Sources/Libraries
################################################################################
output_dir = paste0("gp_experiments/simulations_final/simulations_gpvs/simulated_data")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgpvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgpvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/boxhill_gp.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

# library(fields) # for quilt.plot() ### don't import, it causes problems.
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
xmin = 0
xmax = 1

numx = 21
x_seq = seq(from = xmin, to = xmax, length.out = numx)
x_grid = expand.grid(x_seq, x_seq)

sigmasq_measuremt = 1e-10
sigmasq_signal = 1

################################################################################
# generate functions
registerDoRNG(rng.seed)

if(dimT == 2){
  null_cov = getCov(x_grid, x_grid, typeT, lT, pT, sigmasq_signal) ######################## check
  null_mean = rep(0, numx^2)
} else if(dimT == 1){
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
  output_dir,
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
