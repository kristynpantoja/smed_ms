################################################################################
# last updated: 05/25/2021
# purpose: to test seqmedgp for scenario 3:
#   squared exponential vs. another squared exponential,
#   where the true function is matern

scenario = 4.1 # scenarios: 3.1, 4.1, 5.1, 6.1
input.type = 2 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5

################################################################################
# Sources/Libraries
################################################################################
scenario_subtypes = unlist(strsplit(as.character(scenario), split = "\\."))
output_home = paste0( # works for scenarios1 OR scenarios2 simulations
  "gp_experiments/scenarios", 
  scenario_subtypes[2],
  "/scenario", scenario, "/outputs")
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
sigmasq_measuremt = 1e-10

# shared settings
sigmasq = 1
if(scenario_subtypes[1] %in% c(3, 4)){
  nuggets = c(1e-15, sigmasq_measuremt)
} else if(scenario_subtypes[1] %in% c(5, 6)){
  nuggets = c(1e-5, sigmasq_measuremt)
}
prior_probs = rep(1 / 2, 2)

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
}

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
# read in the data

boxhills = list()
qs = list()
buffers = list()
randoms = list()
spacefills = list()
seqmed.ms = list()

for(i in 1:3){
  # filename_append.tmp for all methods alike
  filename_append.tmp = paste0(
    filename_append, 
    "_input", i, 
    "_seed", rng.seed,
    ".rds"
  )
  boxhills[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_boxhill", 
    filename_append.tmp))
  buffers[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
    "_buffer", 
    "_seq", seq.type,
    filename_append.tmp))
  qs[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
    "_q",
    "_seq", seq.type,
    filename_append.tmp))
  
  randoms[[i]] = readRDS(paste0(
    output_home, 
    "/scenario", scenario, "_random", 
    filename_append.tmp))
  spacefills[[i]] = readRDS(paste0(
    output_home, 
    "/scenario", scenario, "_spacefilling", 
    filename_append.tmp))
  
  seqmed.ms[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
    "_error", 
    "_seq", seq.type,
    filename_append.tmp))
}

################################################################################
# make plots
################################################################################

# input set
bh.in = boxhills[[input.type]]
q.in = qs[[input.type]]
buf.in = buffers[[input.type]]
ran.in = randoms[[input.type]]
sf.in = spacefills[[input.type]]
m.in = seqmed.ms[[input.type]]
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

# all 6 designs
idx = 1
designs = list(bh.in[[idx]], q.in[[idx]], buf.in[[idx]], m.in[[idx]])
design.names = c(
  "bh", "q", "augdist", "measmt")
design.levels = c("measmt", "augdist", "q", "bh")

x.new.mat = matrix(NA, nrow = Nnew, ncol = length(designs))
for(i in 1:length(designs)){
  x.new.mat[, i] = designs[[i]]$x.new
}

data.gg = data.frame(
  index = as.character(rep(1:Nnew, length(designs))), 
  type = factor(rep(design.names, each = Nnew), levels = design.levels), 
  value = as.vector(x.new.mat)
)
data.gg0 = data.frame(
  type = factor(rep(design.names, each = Nin), levels = design.levels), 
  input = rep(x_input, length(designs))
)
text.gg = dplyr::filter(data.gg, index %in% as.character(1:Nnew))
des.plt = ggplot() + 
  geom_point(data = data.gg0, 
             mapping = aes(x = input, y = type)) +
  geom_point(data = data.gg, 
             mapping = aes(x = value, y = type, color = type), 
             inherit.aes = FALSE) + 
  geom_text(data = text.gg, 
            aes(x = value, y = type, label = index), 
            vjust = -0.65 * as.numeric(paste(text.gg$index)), size = 2) +
  xlim(c(xmin, xmax))
des.plt

ggsave(
  filename = paste0("20210525_scen", scenario, "_in", input.type, "_design.pdf"), 
  plot = des.plt, 
  width = 6, height = 4, units = c("in")
)
