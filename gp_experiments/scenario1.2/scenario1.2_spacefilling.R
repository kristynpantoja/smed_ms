################################################################################
# last updated: 05/24/2021
# purpose: to test seqmedgp for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern

scenario = 1.2

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
Nin = 6
numSeq = 15
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10

# space-filling settings

################################################################################
# input data
################################################################################

# 1. make space-filling design
# space_filling = seq(from = xmin, to = xmax, length.out = Nttl)
space_filling_idx = c(1, 1 + ((numx - 1)/(Nttl - 1)) * 1:((numx - 1) / ((numx - 1)/(Nttl - 1))))
space_filling = x_seq[space_filling_idx]

# input set 1 (extrapolation)
x_in1_idx = space_filling_idx[1:Nin]
x_in1 = x_seq[x_in1_idx]
x_spacefill1_idx = space_filling_idx[-c(1:Nin)]
x_spacefill1 = x_seq[x_spacefill1_idx]
# all.equal(space_filling, c(x_in1, x_spacefill1))

# input set 2 (increasing spread)
x_in2_idx = space_filling_idx[c(1, 2, 4, 7, 12, 21)]
x_in2 = x_seq[x_in2_idx]
x_spacefill2_idx = space_filling_idx[-c(1, 2, 4, 7, 12, 21)]
x_spacefill2 = x_seq[x_spacefill2_idx]
# all.equal(space_filling, sort(c(x_in2, x_spacefill2)))

# input set 3 (space-filling / even coverage)
x_in3_idx = c(1, 1 + ((numx - 1)/(Nin - 1)) * 1:((numx - 1) / ((numx - 1)/(Nin - 1))))
x_in3 = x_seq[x_in3_idx]
x_spacefill3_idx = space_filling_idx[!(space_filling_idx %in% x_in3_idx)]
x_spacefill3 = x_seq[x_spacefill3_idx]
# all.equal(space_filling, sort(c(x_in3, x_spacefill3)))

# input set 4 (uniform / random)

################################################################################
# Scenario 1: Squared exponential vs. matern, true = matern
################################################################################
type01 = c("squaredexponential", "matern")
typeT = type01[2]
l01= c(0.01, 0.01)
lT = l01[2]

################################################################################
# import matern functions
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_measuremt), "-")[[1]][2])
}
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario", scenario, "_simulated_functions", filename_append,
  "_seed", rng.seed,
  ".rds"))
numSims = simulated.functions$numSims
x_seq = simulated.functions$x
numx = length(x_seq)
null_cov = simulated.functions$null_cov
null_mean = simulated.functions$null_mean
y_seq_mat = simulated.functions$function_values_mat

################################################################################
# space-filling design

for(j in 1:3){
  
  # j : input setting
  input.type = j
  # input set
  if(input.type == 1){
    x_input = x_in1
    x_input_idx = x_in1_idx
    # new points
    x.new = x_spacefill1
    x.new.idx = x_spacefill1_idx
  } else if(input.type == 2){
    x_input = x_in2
    x_input_idx = x_in2_idx
    # new points
    x.new = x_spacefill2
    x.new.idx = x_spacefill2_idx
  } else if(input.type == 3){
    x_input = x_in3
    x_input_idx = x_in3_idx
    # new points
    x.new = x_spacefill3
    x.new.idx = x_spacefill3_idx
  }
  
  # simulations!
  registerDoRNG(rng.seed)
  spacefills = foreach(
    i = 1:numSims
  ) %dorng% {
    y_seq = y_seq_mat[ , i]
    y_input = y_seq[x_input_idx]
    # new points' y
    y.new = y_seq[x.new.idx]
    list(x = x_input, x.idx = x_input_idx, y = y_input, 
         x.new = x.new, x.new.idx = x.new.idx, y.new = y.new, 
         function.values = y_seq)
  }
  
  filename_append.tmp = filename_append
  filename_append.tmp = paste0(
    filename_append.tmp, 
    "_input", input.type, 
    "_seed", rng.seed,
    ".rds"
  )
  saveRDS(spacefills, 
          file = paste0(
            output_home,
            "/scenario", scenario, "_spacefilling", 
            filename_append.tmp))
}
