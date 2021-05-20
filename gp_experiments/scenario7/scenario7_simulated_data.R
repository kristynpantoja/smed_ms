################################################################################
# last updated: 05/20/2021
# purpose: to test seqmedgp for scenario 7:
#   squared exponential vs. periodic,
#   where the true function is periodic
#   (as opposed to comparing a matern to a periodic, like in scenario 2)

scenario = 7

################################################################################
# Sources/Libraries
################################################################################
output_home = paste0("gp_experiments/scenario", scenario, "/outputs")
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
sigmasq_err = 1e-10

################################################################################
# Scenario 7: Squared exponential vs. periodic, true = periodic
################################################################################
type01 = c("squaredexponential", "periodic")
typeT = type01[2]
l01= c(0.01, 0.01)
lT = l01[2]

################################################################################
# generate periodic functions 
registerDoRNG(rng.seed)
null_cov = getCov(x_seq, x_seq, typeT, lT)
null_mean = rep(0, numx)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values
filename_append = ""
if(!is.null(sigmasq_err)){
  y_seq_mat = y_seq_mat + matrix(rnorm(numx * numSims), 
                                 nrow = nrow(y_seq_mat), ncol = ncol(y_seq_mat))
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_err), "-")[[1]][2])
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
    "/scenario", scenario, "_simulated_functions", filename_append, 
    "_seed", rng.seed,
    ".rds"))
