################################################################################
# last updated: 05/03/2021
# purpose: to test seqmedgp for scenario 4:
#   matern vs. squared exponential,
#   where the true function is periodic

################################################################################
# Sources/Libraries
################################################################################
output_home = "gp_experiments/seqmed_scenario4_20210503/outputs"
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
# errorvar.type = 1 # 1 = phi0 with nugget, 2 = phi1 with nugget
# signalvar.type = 2 # 1 = phi0 sigmasq != 1, 2 = phi1 sigmasq != 1
# input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5

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

# SeqMED settings
sigmasqs = c(1 - 1e-10, 1)
nuggets = c(1e-10, 1e-15)
nugget.sm = NULL
buffer = 0

# boxhill settings
nugget.bh = NULL #1e-10
prior_probs = rep(1 / 2, 2)

# shared settings
sigmasq = 1

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
# Scenario 4: Matern vs. squared exponential, true = periodic
################################################################################
type01 = c("matern", "squaredexponential")
typeT = "periodic"
l01= c(0.01, 0.01)
lT = 0.01

################################################################################
# models - BoxHill
model0.bh = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                 error.var = nugget.bh)
model1.bh = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                 error.var = nugget.bh)

# models - q, random, space-filling, buffer
model0.other = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                    error.var = nugget.sm)
model1.other = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                    error.var = nugget.sm)

# models - SeqMED with different nugget term
# errorvar.type == 1
model0.n1 = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                 error.var = nuggets[1])
model1.n1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                 error.var = nuggets[2])
# errorvar.type == 2
model0.n2 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
                 error.var = nuggets[2])
model1.n2 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                 error.var = nuggets[1])

# models - SeqMED with different signal variance
# signalvar.type == 1
model0.s1 = list(type = type01[1], l = l01[1], signal.var = sigmasqs[1], 
                 error.var = nugget.sm)
model1.s1 = list(type = type01[2], l = l01[2], signal.var = sigmasqs[2], 
                 error.var = nugget.sm)
# signalvar.type == 2
model0.s2 = list(type = type01[1], l = l01[1], signal.var = sigmasqs[2],
                 error.var = nugget.sm)
model1.s2 = list(type = type01[2], l = l01[2], signal.var = sigmasqs[1], 
                 error.var = nugget.sm)

################################################################################
# import periodic functions
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario4_simulated_functions", 
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

boxhills1 = readRDS(paste0(
  output_home, 
  "/scenario4_boxhill", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
boxhills2 = readRDS(paste0(
  output_home, 
  "/scenario4_boxhill", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
boxhills3 = readRDS(paste0(
  output_home, 
  "/scenario4_boxhill", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

qs1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed",
  "_obj", 2,
  "_input", 1,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
qs2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed",
  "_obj", 2,
  "_input", 2,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
qs3 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed",
  "_obj", 2,
  "_input", 3,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

buffers1 = readRDS(paste0(
  output_home,
  "/scenario4_buffer",
  "_obj", 1,
  "_input", 1,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
buffers2 = readRDS(paste0(
  output_home,
  "/scenario4_buffer",
  "_obj", 1,
  "_input", 2,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
buffers3 = readRDS(paste0(
  output_home,
  "/scenario4_buffer",
  "_obj", 1,
  "_input", 3,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

randoms1 = readRDS(paste0(
  output_home, 
  "/scenario4_random", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
randoms2 = readRDS(paste0(
  output_home, 
  "/scenario4_random", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
randoms3 = readRDS(paste0(
  output_home, 
  "/scenario4_random", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

spacefills1 = readRDS(paste0(
  output_home, 
  "/scenario4_spacefilling", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
spacefills2 = readRDS(paste0(
  output_home, 
  "/scenario4_spacefilling", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
spacefills3 = readRDS(paste0(
  output_home, 
  "/scenario4_spacefilling", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

seqmeds.n1.1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n1.2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n1.3 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n2.1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n2.2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n2.3 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s1.1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s1.2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s1.3 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s2.1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s2.2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s2.3 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

################################################################################
# make plots
################################################################################

# models
modelT = list(type = typeT, l = lT, signal.var = sigmasq, error.var = NULL)

boxhills = list(boxhills1, boxhills2, boxhills3)
qs = list(qs1, qs2, qs3)
buffers = list(buffers1, buffers2, buffers3)
randoms = list(randoms1, randoms2, randoms3)
spacefills = list(spacefills1, spacefills2, spacefills3)
seqmed.n1s = list(seqmeds.n1.1, seqmeds.n1.2, seqmeds.n1.3)
seqmed.n2s = list(seqmeds.n2.1, seqmeds.n2.2, seqmeds.n2.3)
seqmed.s1s = list(seqmeds.s1.1, seqmeds.s1.2, seqmeds.s1.3)
seqmed.s2s = list(seqmeds.s2.1, seqmeds.s2.2, seqmeds.s2.3)

# calculate the RSS01
getRSST = function(design, model){
  pred.tmp = getGPPredictive(design$x.new, design$x, design$y, 
                             model$type, model$l, 
                             model$signal.var, model$error.var)
  RSST.tmp = sum((pred.tmp$pred_mean - design$y.new)^2, na.rm = TRUE) 
  return(data.frame("RSST" = RSST.tmp))
}

RSS.df = data.frame(RSST = numeric(), 
                    type = character(), sim = numeric(), input = numeric())
for(k in 1:3){
  for(j in 1:numSims){
    # designs at sim b
    bh = boxhills[[k]][[j]]
    q = qs[[k]][[j]]
    b = buffers[[k]][[j]]
    r = randoms[[k]][[j]]
    sf = spacefills[[k]][[j]]
    n1 = seqmed.n1s[[k]][[j]]
    n2 = seqmed.n2s[[k]][[j]]
    s1 = seqmed.s1s[[k]][[j]]
    s2 = seqmed.s2s[[k]][[j]]
    # sequence of PPHs for each design
    RSS.bh = getRSST(bh, modelT)
    RSS.q = getRSST(q, modelT)
    RSS.b = getRSST(b, modelT)
    RSS.r = getRSST(r, modelT)
    RSS.sf = getRSST(sf, modelT)
    RSS.n1 = getRSST(n1, modelT)
    RSS.n2 = getRSST(n2, modelT)
    RSS.s1 = getRSST(s1, modelT)
    RSS.s2 = getRSST(s2, modelT)
    # master data frame
    RSS.bh$type = "boxhill"
    RSS.q$type = "q"
    RSS.b$type = "buffer"
    RSS.r$type = "random"
    RSS.sf$type = "spacefill"
    RSS.n1$type = "nugget1"
    RSS.n2$type = "nugget2"
    RSS.s1$type = "signal1"
    RSS.s2$type = "signal2"
    RSS.tmp = rbind(
      RSS.bh, RSS.q, RSS.b, RSS.r, RSS.sf, 
      RSS.n1, RSS.n2, RSS.s1, RSS.s2)
    RSS.tmp$sim = j
    RSS.tmp$input = k
    RSS.df = rbind(RSS.df, RSS.tmp)
  }
}

RSSTmean = aggregate(RSS.df$RSST, by = list(RSS.df$type, RSS.df$input), 
                     FUN = function(x) mean(x, na.rm = TRUE))
names(RSSTmean) = c("type", "input", "value")

RSSTmean$type = factor(RSSTmean$type)
RSSTmean$input = factor(RSSTmean$input, 
                        labels = c("extrapolation", "inc spread", "even coverage"))

# RSST
RSST.plt = ggplot(RSSTmean, 
                  aes(x = input, y = value, group = type, color = type)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "RSST", x = "Initial Data")
RSST.plt
