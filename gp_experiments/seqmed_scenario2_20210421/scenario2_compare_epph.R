################################################################################
# last updated: 04/21/2021
# purpose: to test seqmedgp for scenario 2:
#   matern vs. periodic,
#   where the true function is periodic

################################################################################
# Sources/Libraries
################################################################################
output_home = "gp_experiments/seqmed_scenario2_20210421/outputs"
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
input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
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
nuggets = c(1e-5, 1e-10) # had to change, to fix solve() issue
nugget.q = NULL # nugget for q, random, and space-filling designs
nugget.sm = 1e-10 # nugget for different signal variances and buffer
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
# Scenario 2: Matern vs. periodic, true = periodic
################################################################################
type01 = c("matern", "periodic")
typeT = type01[2]
l01= c(0.01, 0.1)
lT = l01[2]

################################################################################
# models - BoxHill
model0.bh = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                 error.var = nugget.bh)
model1.bh = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                 error.var = nugget.bh)

# models - q, random, space-filling, \sout{buffer}
model0.q = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                    error.var = nugget.q)
model1.q = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                    error.var = nugget.q)

# models - different buffer
model0.sm = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                    error.var = nugget.sm)
model1.sm = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
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
  "/scenario2_simulated_functions", 
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

boxhills = readRDS(paste0(
  output_home, 
  "/scenario2_boxhill", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

qs = readRDS(paste0(
  output_home,
  "/scenario2_seqmed",
  "_obj", 2,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

buffers = readRDS(paste0(
  output_home,
  "/scenario2_buffer",
  "_obj", 1,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

randoms = readRDS(paste0(
  output_home, 
  "/scenario2_random", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

spacefills = readRDS(paste0(
  output_home, 
  "/scenario2_spacefilling", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

seqmeds.n1 = readRDS(paste0(
  output_home,
  "/scenario2_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n2 = readRDS(paste0(
  output_home,
  "/scenario2_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s1 = readRDS(paste0(
  output_home,
  "/scenario2_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s2 = readRDS(paste0(
  output_home,
  "/scenario2_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

################################################################################
# make plots
################################################################################
PPHs_seq = list()

# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
              error.var = NULL)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = NULL)

getPPHseq = function(design, model0, model1){
  PPH0_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  PPH1_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  for(i in 1:length(as.vector(na.omit(design$y.new)))){
    y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:i])
    x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:i])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1)
      )
    )
    PPH0_seq[i] = PPHs.tmp[1]
    PPH1_seq[i] = PPHs.tmp[2]
  }
  if(length(PPH0_seq) < Nnew){
    PPH0_seq[(length(PPH0_seq) + 1):Nnew] = NA
    PPH1_seq[(length(PPH1_seq) + 1):Nnew] = NA
  }
  return(data.frame(
    index = 1:Nnew, 
    PPH0 = PPH0_seq, 
    PPH1 = PPH1_seq
  ))
}

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  type = character(), sim = numeric())
for(j in 1:numSims){
  # designs at sim b
  bh = boxhills[[j]]
  q = qs[[j]]
  b = buffers[[j]]
  r = randoms[[j]]
  sf = spacefills[[j]]
  n1 = seqmeds.n1[[j]]
  n2 = seqmeds.n2[[j]]
  s1 = seqmeds.s1[[j]]
  s2 = seqmeds.s2[[j]]
  # sequence of PPHs for each design
  PPH_seq.bh = getPPHseq(bh, model0, model1) #model0.bh, model1.bh)
  PPH_seq.q = getPPHseq(q, model0, model1) #model0.q, model1.q)
  PPH_seq.b = getPPHseq(b, model0, model1) #model0.sm, model1.sm)
  PPH_seq.r = getPPHseq(r, model0, model1) #model0.q, model1.q)
  PPH_seq.sf = getPPHseq(sf, model0, model1) #model0.q, model1.q)
  PPH_seq.n1 = getPPHseq(n1, model0, model1) #model0.n1, model1.n1)
  PPH_seq.n2 = getPPHseq(n2, model0, model1) #model0.n2, model1.n2)
  PPH_seq.s1 = getPPHseq(s1, model0, model1) #model0.s1, model1.s1)
  PPH_seq.s2 = getPPHseq(s2, model0, model1) #model0.s2, model1.s2)
  # master data frame
  PPH_seq.bh$type = "boxhill"
  PPH_seq.q$type = "q"
  PPH_seq.b$type = "buffer"
  PPH_seq.r$type = "random"
  PPH_seq.sf$type = "spacefill"
  PPH_seq.n1$type = "nugget1"
  PPH_seq.n2$type = "nugget2"
  PPH_seq.s1$type = "signal1"
  PPH_seq.s2$type = "signal2"
  PPH_seq.tmp = rbind(
    PPH_seq.bh, PPH_seq.q, PPH_seq.b, PPH_seq.r, PPH_seq.sf, 
    PPH_seq.n1, PPH_seq.n2, PPH_seq.s1, PPH_seq.s2)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPH0mean_seq = aggregate(PPH_seq$PPH0, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean_seq) = c("index", "type", "value")
PPH0mean_seq$Hypothesis = "H0"
PPH1mean_seq = aggregate(PPH_seq$PPH1, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean_seq) = c("index", "type", "value")
PPH1mean_seq$Hypothesis = "H1"

PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq)
ggplot(PPHmean_seq, aes(x = index, y = value, color = type, linetype = type)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  theme_bw() +
  ylim(0, 1)
