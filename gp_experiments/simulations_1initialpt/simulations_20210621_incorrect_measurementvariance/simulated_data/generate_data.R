################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 1 or 2
#   where H1 is true

typeT = "periodic"
lT = 0.01

################################################################################
# Sources/Libraries
################################################################################
output_home = paste0("gp_experiments/simulated_data")
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMEDgp.R", sep = ""))
source(paste(functions_home, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_home, "/boxhill.R", sep = ""))
source(paste(functions_home, "/boxhill_gp.R", sep = ""))
source(paste(functions_home, "/kl_divergence.R", sep = ""))

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
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

################################################################################
# generate functions 
registerDoRNG(rng.seed)
null_cov = getCov(x_seq, x_seq, typeT, lT, sigmasq_signal)
null_mean = rep(0, numx)

# the function values
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) 
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  y_seq_mat = y_seq_mat + 
    matrix(rnorm(numx * numSims, 1, sqrt(sigmasq_measuremt)), 
           nrow = nrow(y_seq_mat), ncol = ncol(y_seq_mat))
  filename_append = "_noise"
    # paste0("_noise", strsplit(as.character(sigmasq_measuremt), "-")[[1]][2])
}

saveRDS(
  list(
    x = x_seq, 
    null_mean = null_mean, 
    null_cov = null_cov, 
    numSims = numSims, 
    function_values_mat = y_seq_mat
  ), 
  file = paste0(
    output_home,
    "/", typeT,
    "_l", lT,
    filename_append, 
    "_seed", rng.seed,
    ".rds"))
