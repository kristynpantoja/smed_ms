################################################################################
# last updated: 04/21/2021
# purpose: to test seqmedgp for scenario 2:
#   matern vs. periodic,
#   where the true function is periodic

################################################################################
# Sources/Libraries
################################################################################
output_home = "gp_experiments/seqmed_scenario2_20210421/outputs"
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

################################################################################
# Scenario 2: Matern vs. periodic, true = periodic
################################################################################
type01 = c("matern", "periodic")
typeT = type01[2]
l01= c(0.1, 0.5)
lT = l01[2]

################################################################################
# generate periodic functions 
registerDoRNG(rng.seed)
null_cov = getCov(x_seq, x_seq, typeT, lT)
null_mean = rep(0, numx)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values

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
    "/scenario2_simulated_functions", 
    "_seed", rng.seed,
    ".rds"))
