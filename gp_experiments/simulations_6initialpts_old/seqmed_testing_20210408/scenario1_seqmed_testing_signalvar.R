################################################################################
# last updated: 04/14/2021
# purpose: to test seqmedgp for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern
# trying out some (not necessarily MED) designs
# changed SeqMEDgp to take in model0, model1

################################################################################
# Sources/Libraries
################################################################################
output_home = "run_designs/gp_experiments"
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
signalvar.type = 1 # 1 = phi0 sigmasq != 1, 2 = phi1 sigmasq != 1
input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5

# simulations settings
numSims = 10
Nin = 6
if(seq.type == 1){
  numSeq = 15
  seqN = 1
} else if(seq.type == 2){
  numSeq = 3
  seqN = 5
}
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)

# SeqMED settings
sigmasqs = c(1 - 1e-15, 1)
nuggetSM = NULL
buffer = 0

# boxhill settings
prior_probs = rep(1 / 2, 2)
nuggetBH = 1e-10

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
l01= c(0.01, 0.01) # SIM SETTING
# l01= c(0.1, 0.1) # DEMO SETTING

# generate matern functions ####################################################
registerDoRNG(rng.seed)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values

# bh settings
if(signalvar.type == 1){
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasqs[1], 
                error.var = nuggetSM)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasqs[2], 
                error.var = nuggetSM)
  
} else if(signalvar.type == 2){
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasqs[2],
                error.var = nuggetSM)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasqs[1], 
                error.var = nuggetSM)
  
}

################################################################################
# generate seqmeds #############################################################
################################################################################

################################################################################
# try different objective.type - SIM SETTINGS ##################################
################################################################################

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

# seqmed
registerDoRNG(rng.seed)
seqmeds = foreach(
  i = 1:numSims
) %dorng% {
  y_seq = y_seq_mat[ , i]
  y_input = y_seq[x_input_idx]
  SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, model0 = model0, model1 = model1, 
    numSeq = numSeq, seqN = seqN, prints = TRUE, buffer = buffer, 
    objective.type = 1, seed = 1234)
}