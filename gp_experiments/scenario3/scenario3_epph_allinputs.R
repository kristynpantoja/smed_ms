################################################################################
# last updated: 04/28/2021
# purpose: to test seqmedgp for scenario 3:
#   squared exponential vs. another squared exponential,
#   where the true function is matern

################################################################################
# Sources/Libraries
################################################################################
output_home = "gp_experiments/seqmed_scenario3_20210416/outputs"
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
sigmasq_err = 1e-10

# SeqMED settings
sigmasqs = c(1 - 1e-10, 1)
nuggets = c(1e-5, 1e-10)
buffer = 0

# boxhill settings
prior_probs = rep(1 / 2, 2)

# shared settings
sigmasq = 1
nugget = sigmasq_err

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
# Scenario 3: Squared exponential vs. squared exponential, true = matern
################################################################################
type01 = c("squaredexponential", "squaredexponential")
typeT = "matern"
l01= c(0.005, 0.01)
lT = 0.01

################################################################################
# import matern functions
filename_append = ""
if(!is.null(sigmasq_err)){
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_err), "-")[[1]][2])
}
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario3_simulated_functions", filename_append,
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
  "/scenario3_boxhill", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
boxhills2 = readRDS(paste0(
  output_home, 
  "/scenario3_boxhill", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
boxhills3 = readRDS(paste0(
  output_home, 
  "/scenario3_boxhill", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

qs1 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed",
  "_obj", 2,
  "_input", 1,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
qs2 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed",
  "_obj", 2,
  "_input", 2,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
qs3 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed",
  "_obj", 2,
  "_input", 3,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

buffers1 = readRDS(paste0(
  output_home,
  "/scenario3_buffer",
  "_obj", 1,
  "_input", 1,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
buffers2 = readRDS(paste0(
  output_home,
  "/scenario3_buffer",
  "_obj", 1,
  "_input", 2,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
buffers3 = readRDS(paste0(
  output_home,
  "/scenario3_buffer",
  "_obj", 1,
  "_input", 3,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

randoms1 = readRDS(paste0(
  output_home, 
  "/scenario3_random", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
randoms2 = readRDS(paste0(
  output_home, 
  "/scenario3_random", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
randoms3 = readRDS(paste0(
  output_home, 
  "/scenario3_random", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

spacefills1 = readRDS(paste0(
  output_home, 
  "/scenario3_spacefilling", 
  "_input", 1, 
  "_seed", rng.seed, 
  ".rds"
))
spacefills2 = readRDS(paste0(
  output_home, 
  "/scenario3_spacefilling", 
  "_input", 2, 
  "_seed", rng.seed, 
  ".rds"
))
spacefills3 = readRDS(paste0(
  output_home, 
  "/scenario3_spacefilling", 
  "_input", 3, 
  "_seed", rng.seed, 
  ".rds"
))

seqmeds.n1.1 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n1.2 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n1.3 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n2.1 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n2.2 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.n2.3 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s1.1 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s1.2 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s1.3 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", 3, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s2.1 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", 1, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s2.2 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", 2, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))
seqmeds.s2.3 = readRDS(paste0(
  output_home,
  "/scenario3_seqmed", 
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
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
              error.var = NULL)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = NULL)
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

# calculate the final posterior probability
getPPH = function(design, model0, model1, modelT){
  y.tmp = c(design$y, as.vector(na.omit(design$y.new)))
  x.tmp = c(design$x, as.vector(na.omit(design$x.new)))
  PPHs.tmp = getHypothesesPosteriors(
    prior.probs = rep(1 / 3, 3), 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0),
      Evidence_gp(y.tmp, x.tmp, model1),
      Evidence_gp(y.tmp, x.tmp, modelT)
    )
  )
  return(data.frame("H0" = PPHs.tmp[1], "H1" = PPHs.tmp[2], "HT" = PPHs.tmp[3]))
}

# true model
modelT = list(type = typeT, l = lT, signal.var = sigmasq, error.var = NULL)

PPH = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
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
    PPH.bh = getPPH(bh, model0, model1, modelT) # model0.bh, model1.bh, modelT)
    PPH.q = getPPH(q, model0, model1, modelT) # model0.other, model1.other, modelT)
    PPH.b = getPPH(b, model0, model1, modelT) # model0.other, model1.other, modelT)
    PPH.r = getPPH(r, model0, model1, modelT) # model0.other, model1.other, modelT)
    PPH.sf = getPPH(sf, model0, model1, modelT) # model0.other, model1.other, modelT)
    PPH.n1 = getPPH(n1, model0, model1, modelT) # model0.n1, model1.n1, modelT)
    PPH.n2 = getPPH(n2, model0, model1, modelT) # model0.n2, model1.n2, modelT)
    PPH.s1 = getPPH(s1, model0, model1, modelT) # model0.s1, model1.s1, modelT)
    PPH.s2 = getPPH(s2, model0, model1, modelT) # model0.s2, model1.s2, modelT)
    # master data frame
    PPH.bh$type = "boxhill"
    PPH.q$type = "q"
    PPH.b$type = "buffer"
    PPH.r$type = "random"
    PPH.sf$type = "spacefill"
    PPH.n1$type = "nugget1"
    PPH.n2$type = "nugget2"
    PPH.s1$type = "signal1"
    PPH.s2$type = "signal2"
    PPH.tmp = rbind(
      PPH.bh, PPH.q, PPH.b, PPH.r, PPH.sf, 
      PPH.n1, PPH.n2, PPH.s1, PPH.s2)
    PPH.tmp$sim = j
    PPH.tmp$input = k
    PPH = rbind(PPH, PPH.tmp)
  }
}

PPH0mean = aggregate(PPH$H0, by = list(PPH$type, PPH$input), 
                     FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean) = c("type", "input", "value")
PPH0mean$Hypothesis = "H0"
PPH1mean = aggregate(PPH$H1, by = list(PPH$type, PPH$input), 
                     FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean) = c("type", "input", "value")
PPH1mean$Hypothesis = "H1"
PPHTmean = aggregate(PPH$HT, by = list(PPH$type, PPH$input), 
                     FUN = function(x) mean(x, na.rm = TRUE))
names(PPHTmean) = c("type", "input", "value")
PPHTmean$Hypothesis = "HT"

PPHmean = rbind(PPH0mean, PPH1mean, PPHTmean)
PPHmean$type = factor(PPHmean$type)
PPHmean$Hypothesis = factor(PPHmean$Hypothesis)
PPHmean$input = factor(PPHmean$input, 
                       labels = c("extrapolation", "inc spread", "even coverage"))

PPHT.plt = ggplot(dplyr::filter(PPHmean, Hypothesis == "HT"), 
                  aes(x = input, y = value, group = type, color = type)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  ylim(0, 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "P(HT|X, Y)", x = "Initial Data")
PPHT.plt
