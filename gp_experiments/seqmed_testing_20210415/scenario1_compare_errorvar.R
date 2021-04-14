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
output_home = "gp_experiments/seqmed_testing_20210415/outputs"
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
nuggets = c(1e-10, 1e-15)
buffer = 0

# boxhill settings
nugget.bh = 1e-10
prior_probs = rep(1 / 2, 2)

# shared settings
sigmasq = 1

# models - BoxHill
model0.bh = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
              error.var = nugget)
model1.bh = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = nugget)

# models - SeqMED
# errorvar.type == 1
model0.n1 = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
              error.var = nuggets[1])
model1.n1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = nuggets[2])
#errorvar.type == 2
model0.n2 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
              error.var = nuggets[2])
model1.n2 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = nuggets[1])

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

################################################################################
# read in the data
################################################################################

boxhill = readRDS(paste0(
  output_home, 
  "/scenario1_boxhill", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

seqmed1 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
  "_nugget", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmed2 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
  "_nugget", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

################################################################################
# posterior probabilities
################################################################################

PPHs = list()
PPHs_seq = list()

idx = 1
bh = boxhill[[2]][[idx]]
s1 = seqmed1[[2]][[idx]]
s2 = seqmed2[[2]][[idx]]

# PPHs for BH
y.tmp = c(bh$y, as.vector(na.omit(bh$y.new)))
x.tmp = c(bh$x, as.vector(na.omit(bh$x.new)))
PPHs[[1]] = getHypothesesPosteriors(
  prior.probs = prior_probs, 
  evidences = c(
    Evidence_gp(y.tmp, x.tmp, model0),
    Evidence_gp(y.tmp, x.tmp, model1)
  )
)

# PPHs_seq for BH
PPH0_seq.tmp = rep(NA, length(as.vector(na.omit(bh$y.new))))
PPH1_seq.tmp = rep(NA, length(as.vector(na.omit(bh$y.new))))
for(i in 1:length(as.vector(na.omit(bh$y.new)))){
  y.tmp = c(bh$y, as.vector(na.omit(bh$y.new))[1:i])
  x.tmp = c(bh$x, as.vector(na.omit(bh$x.new))[1:i])
  PPHs.tmp = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )
  PPH0_seq.tmp[i] = PPHs.tmp[1]
  PPH1_seq.tmp[i] = PPHs.tmp[2]
}
if(length(PPH0_seq.tmp) < Nnew){
  PPH0_seq.tmp[(length(PPH0_seq.tmp) + 1):Nnew] = NA
  PPH1_seq.tmp[(length(PPH1_seq.tmp) + 1):Nnew] = NA
}
PPHs_seq[[1]] = data.frame(
  "value" = c(PPH0_seq.tmp, PPH1_seq.tmp), 
  "Hypothesis" = c(rep("H0", Nnew), rep("H1", Nnew))
)

# PPHs for SeqMEDs
for(k in 1:length(seqmed_list)){
  # PPHs for SeqMED[[k]]
  y.tmp = c(seqmed_list[[k]]$y, seqmed_list[[k]]$y.new)
  x.tmp = c(seqmed_list[[k]]$x, seqmed_list[[k]]$x.new)
  PPHs[[k + 1]] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
        Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
      )
    )
  )
  # PPHs_seq for SeqMED[[k]]
  PPH0_seq.tmp = rep(NA, Nnew)
  PPH1_seq.tmp = rep(NA, Nnew)
  for(i in 1:Nnew){
    y.tmp2 = c(seqmed_list[[k]]$y, seqmed_list[[k]]$y.new[1:i])
    x.tmp2 = c(seqmed_list[[k]]$x, seqmed_list[[k]]$x.new[1:i])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gp(y.tmp2, x.tmp2, model0, nuggetSM),
        Evidence_gp(y.tmp2, x.tmp2, model1, nuggetSM)
      )
    )
    PPH0_seq.tmp[i] = PPHs.tmp[1]
    PPH1_seq.tmp[i] = PPHs.tmp[2]
  }
  # if(length(PPH0_seq.tmp) < Nnew){
  #   PPH0_seq.tmp[(length(PPH0_seq.tmp) + 1):Nnew] = NA
  #   PPH1_seq.tmp[(length(PPH1_seq.tmp) + 1):Nnew] = NA
  # }
  PPHs_seq[[k + 1]] = data.frame(
    "value" = c(PPH0_seq.tmp, PPH1_seq.tmp), 
    "Hypothesis" = c(rep("H0", Nnew), rep("H1", Nnew))
  )
}

# melt PPHs and PPHs_seq
# PPHs
PPHs.df = matrix(NA, nrow = length(types), ncol = 2)
for(k in 1:length(PPHs)){
  PPHs.df[k, ] = PPHs[[k]]
}
PPHs.df = as.data.frame(PPHs.df)
rownames(PPHs.df) = types
names(PPHs.df) = c("H0", "H1")
PPHs.df = melt(PPHs.df)
PPHs.df$Type = rep(types, 2)
colnames(PPHs.df)[1] = "Hypothesis"
# PPHs_seq
PPHs_seq.df = PPHs_seq[[1]]
PPHs_seq.df$Type = types[1]
PPHs_seq.df$Index = rep(1:Nnew, 2)
for(k in 2:length(PPHs_seq)){
  PPHs_seq_k.tmp = PPHs_seq[[k]]
  PPHs_seq_k.tmp$Type = types[k]
  PPHs_seq_k.tmp$Index = rep(1:Nnew, 2)
  PPHs_seq.df = rbind(PPHs_seq.df, PPHs_seq_k.tmp)
}

# plot PPHs
ggplot(PPHs.df, aes(x = Type, y = value, color = Type)) + 
  facet_wrap(~Hypothesis) + 
  geom_point()

# plot PPHs_seq
ggplot(PPHs_seq.df, aes(x = Index, y = value, color = Type)) + 
  facet_wrap(~Hypothesis) + 
  geom_path()

# results are very strange. 
#   MED doesn't do well unless it's in stages.



################################################################################
# posterior probabilities for each possible design point #######################
################################################################################

# function for calculating PPH, given initial data
getPPH = function(x.idx, which.pph = 0){
  y.tmp = c(bh$y, y_seq[x.idx])
  x.tmp = c(bh$x, x_seq[x.idx])
  PPHs = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )
  if(which.pph == 0){
    return(PPHs[1])
  } else if(which.pph == 1){
    return(PPHs[2])
  } else{
    stop("getPPH: invalid which.pph argument -- must be 0 or 1")
  }
}

x_idx_seq = 1:numx
pph0_seq = sapply(x_idx_seq, FUN = function(idx) getPPH(idx, 0))
pph1_seq = sapply(x_idx_seq, FUN = function(idx) getPPH(idx, 1))

pph0_seq.df = data.frame(
  Hypothesis = "H0", value = pph0_seq, x = x_seq, Index = x_idx_seq)
pph1_seq.df = data.frame(
  Hypothesis = "H1", value = pph1_seq, x = x_seq, Index = x_idx_seq)
pph_seq.df = rbind(pph0_seq.df, pph1_seq.df)
ggplot(pph_seq.df, aes(x = x, y = value, color = Hypothesis)) + 
  geom_vline(xintercept = x_input, color = "gray") +
  # geom_vline(xintercept = x.new.mat[, "q"], color = "lightgray") +
  geom_path() + 
  geom_point(data = data.frame(x = x_input, y = 0.5), 
             mapping = aes(x = x, y = y), 
             inherit.aes = FALSE)
