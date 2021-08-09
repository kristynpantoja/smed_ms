################################################################################
# last updated: 07/13/2021
# purpose: to make random design for all types of data

typeT = "squaredexponential"
lT = 0.01
dimT = 2 #1, 2
# if 1, assume 1st dimension is true.

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_vs"
output_dir = paste0(sims_dir, "/fixed_designs/outputs")
data_dir = paste0(sims_dir, "/simulated_data")
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
# numSims = 25
Nin = 3
numSeq = 9
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = 0
xmax = 1

sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# random settings

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
simulated_data_file = paste0(
  data_dir,
  "/", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append, 
  "_seed", rng.seed,
  ".rds")
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x_seq
x_grid = as.matrix(simulated.data$x_grid)
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = sample(1:numx, Nin)
x_input = x_grid[x_input_idx + numx * (1:Nin), , drop = FALSE]

################################################################################
# generate random designs (sample x ~ U[xmin, xmax])

# simulations!
registerDoRNG(rng.seed)
randoms = foreach(
  i = 1:numSims
) %dorng% {
  y_seq = y_seq_mat[ , i]
  y_input = y_seq[x_input_idx]
  
  x.new.idx  = sample(1:dim(x_grid)[1], Nnew)
  x.new = x_grid[x.new.idx, ]
  y.new = y_seq[x.new.idx]
  list(x = x_input, x.idx = x_input_idx, y = y_input, 
       x.new = x.new, x.new.idx = x.new.idx, y.new = y.new, 
       function.values = y_seq)
}

filename_append.tmp = filename_append
filename_append.tmp = paste0(
  filename_append.tmp, 
  "_seed", rng.seed,
  ".rds"
)
simulated_spacefilling_file = paste0(
  output_dir,
  "/random", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT, 
  filename_append.tmp)
saveRDS(randoms, file = simulated_spacefilling_file)
