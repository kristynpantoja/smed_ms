################################################################################
# last updated: 07/13/2021
# purpose: to test SeqMEDgpvs()

dimT = 2 # 1, 2
seq.type = 1
lT = 0.01
typeT = "squaredexponential"

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/gpvs"
output_dir = paste0(sims_dir, "/modelselection_designs/outputs")
data_dir = paste0(sims_dir, "/simulated_data/outputs")
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
rng.seed = 123

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

# library(expm)
# library(matrixStats)
# library(MASS)
# library(mvtnorm)
# library(knitr)

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
Nin = 1
if(seq.type == 1){
  numSeq = 9
  seqN = 1
} else if(seq.type == 2){
  numSeq = 3
  seqN = 3
}
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = 0
xmax = 1
p = 2
k = 4 * p
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# seqmed settings
obj_type = 4
newq = FALSE

################################################################################
# Scenario settings
################################################################################

l01= c(lT, lT)
type01 = c(typeT, typeT)
indices0 = c(1)
indices1 = c(1, 2)

model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
              indices = indices0,
              measurement.var = sigmasq_measuremt)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
              indices = indices1, 
              measurement.var = sigmasq_measuremt)

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
# generate designs

# simulations!
registerDoRNG(rng.seed)
seqmeds = foreach(
  i = 1:numSims
) %dorng% {
  y_seq = y_seq_mat[ , i]
  y_input = y_seq[x_input_idx]
  
  SeqMEDgpvs(
    y.in = y_input, x.in = x_input, x.in.idx = x_input_idx,
    candidates = x_grid, function.values = y_seq, 
    xmin = xmin, xmax = xmax, k = k, p = p, 
    numSeq = numSeq, seqN = seqN, 
    model0 = model0, model1 = model1, 
    objective.type = obj_type, newq = newq)
}

filename_append.tmp = paste0(
  filename_append, 
  "_seed", rng.seed,
  ".rds"
)
saveRDS(seqmeds, 
        file = paste0(
          output_dir,
          "/seqmed", 
          "_keepq_cap", 
          "_seq", seq.type,
          "_l", lT,
          "_dim", dimT,
          filename_append.tmp))

