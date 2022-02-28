################################################################################
# last updated: 05/27/2021
# purpose: to make random design for all types of data

typeT = "periodic"
pT = 0.05
lT = 0.5
sigmasq_signal_seq = seq(0.5, 1.5, length.out = 101)

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_MSP/scenario_MSP_sigvargrid"
output_dir = paste0(sims_dir, "/spacefilling_designs/outputs")
data_dir = paste0(sims_dir, "/simulated_data")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgp.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgp_batch.R", sep = ""))
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
numSims = 25
Nin = 1
numSeq = 15
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# random settings

################################################################################

# simulations!
for(j in 1:length(sigmasq_signal_seq)){
  sigmasq_signal = sigmasq_signal_seq[j]
  
  # import data
  filename_append = ""
  if(!is.null(sigmasq_measuremt)){
    filename_append = "_noise"
  }
  if(typeT == "periodic"){
    if(is.null(pT)) pT = 0.26
    simulated_data_file = paste0(
      data_dir,
      "/", typeT,
      "_l", lT,
      "_p", pT,
      "_sig", sigmasq_signal,
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  } else{
    simulated_data_file = paste0(
      data_dir,
      "/", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  }
  simulated.data = readRDS(simulated_data_file)
  numSims = simulated.data$numSims
  x_seq = simulated.data$x
  numx = length(x_seq)
  null_cov = simulated.data$null_cov
  null_mean = simulated.data$null_mean
  y_seq_mat = simulated.data$function_values_mat
  
  ################################################################################
  # initial design
  x_input_idx = ceiling(numx / 2)
  x_input = x_seq[x_input_idx]
  
  ################################################################################
  # generate random designs (sample x ~ U[xmin, xmax])
  
  registerDoRNG(rng.seed)
  randoms = foreach(
    i = 1:numSims
  ) %dorng% {
    y_seq = y_seq_mat[ , i]
    y_input = y_seq[x_input_idx]
    
    x.new.idx  = sample(1:numx, Nnew)
    x.new = x_seq[x.new.idx]
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
  if(typeT == "periodic"){
    simulated_spacefilling_file = paste0(
      output_dir,
      "/random", 
      "_", typeT,
      "_l", lT,
      "_p", pT,
      "_sig", sigmasq_signal,
      filename_append.tmp)
  } else{
    simulated_spacefilling_file = paste0(
      output_dir,
      "/random", 
      "_", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
      filename_append.tmp)
  }
  saveRDS(randoms, file = simulated_spacefilling_file)
}