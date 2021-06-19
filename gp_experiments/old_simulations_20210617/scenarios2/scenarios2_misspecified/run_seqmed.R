################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 3.2, 4.2, 5.2, or 6.2 (H0, H1 different sigmasq_signal)
#   where both hypotheses are misspecified

scenario = 6.2
scenario_subtypes = unlist(strsplit(as.character(scenario), split = "\\."))

################################################################################
# Sources/Libraries
################################################################################
output_home = paste0("gp_experiments/scenarios2/scenarios2_misspecified/outputs")
data_home = "gp_experiments/simulated_data"
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
nworkers = detectCores() - 2
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
Nnew = 15
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# SeqMED settings
sigmasq_signals = c(1 - 1e-10, 1)
nugget = sigmasq_measuremt
buffer = 0

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
# Scenario settings
################################################################################
if(scenario_subtypes[1] == 3){
  type01 = c("squaredexponential", "squaredexponential")
  typeT = "matern"
  l01= c(0.005, 0.01)
  lT = 0.01
} else if(scenario_subtypes[1] == 4){
  type01 = c("matern", "squaredexponential")
  typeT = "periodic"
  l01= c(0.01, 0.01)
  lT = 0.01
} else if(scenario_subtypes[1] == 5){
  type01 = c("matern", "periodic")
  typeT = "squaredexponential"
  l01= c(0.01, 0.01)
  lT = 0.01
} else if(scenario_subtypes[1] == 6){
  type01 = c("squaredexponential", "periodic")
  typeT = "matern"
  l01= c(0.01, 0.01)
  lT = 0.01
} else{
  stop("invalid scenario number")
}

################################################################################
# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signals[1], 
              measurement.var = nugget)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signals[2], 
              measurement.var = nugget)

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
simulated.data = readRDS(paste0(
  data_home,
  "/", typeT,
  "_l", lT,
  filename_append, 
  "_seed", rng.seed,
  ".rds"))
numSims = simulated.data$numSims
x_seq = simulated.data$x
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# generate seqmeds 

for(j in 1:3){
  for(k in 1:2){
    
    # j : input setting
    input.type = j
    # input set
    if(input.type == 1){
      x_input = x_in1
      x_input_idx = x_in1_idx
    } else if(input.type == 2){
      x_input = x_in2
      x_input_idx = x_in2_idx
    } else if(input.type == 3){
      x_input = x_in3
      x_input_idx = x_in3_idx
    }
    
    # k : sequential setting
    seq.type = k
    if(seq.type == 1){
      numSeq = 15
      seqN = 1
    } else if(seq.type == 2){
      numSeq = 3
      seqN = 5
    }
    
    # simulations!
    registerDoRNG(rng.seed)
    seqmeds = foreach(
      b = 1:numSims
    ) %dorng% {
      y_seq = y_seq_mat[ , b]
      y_input = y_seq[x_input_idx]
      SeqMEDgp(
        y0 = y_input, x0 = x_input, x0.idx = x_input_idx, 
        candidates = x_seq, function.values = y_seq, 
        model0 = model0, model1 = model1, 
        numSeq = numSeq, seqN = seqN, prints = FALSE, buffer = buffer, 
        objective.type = 1, noise = FALSE, measurement.var = sigmasq_measuremt)
    }
    
    print(paste0("completed j = ", j, ", k = ", k, "!"))
    
    filename_append.tmp = paste0(
      filename_append, 
      "_input", input.type, 
      "_seed", rng.seed,
      ".rds"
    )
    saveRDS(seqmeds, 
            file = paste0(
              output_home,
              "/scenario", scenario, "_seqmed", 
              "_seq", seq.type,
              filename_append.tmp))
  }
}

