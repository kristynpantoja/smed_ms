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

# SeqMED settings
sigmasq = 1
sigmasqs = c(1 - 1e-10, 1)
nuggets = c(1e-10, 1e-15) # had to change, to fix solve() issue
nugget = NULL
nugget.q = nugget # nugget for q, random, and space-filling designs
nugget.sm = nugget # nugget for different signal variances and buffer
buffer = 0

# boxhill settings
nugget.bh = nugget #1e-10
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

boxhills = list()
qs = list()
buffers = list()
randoms = list()
spacefills = list()
seqmed.n1s = list()
seqmed.n2s = list()
seqmed.s1s = list()
seqmed.s2s = list()

for(i in 1:3){
  
  file_name_end0 = paste0(
    "_input", i, 
    "_seed", rng.seed,
    ".rds"
  )
  if(!is.null(nugget)){
    file_name_end = paste0(
      "_nugget", strsplit(as.character(nugget), "-")[[1]][2], 
      file_name_end0)
  } else{
    file_name_end = file_name_end0
  }
  
  boxhills[[i]] = readRDS(paste0(
    output_home, 
    "/scenario4_boxhill", 
    file_name_end))
  qs[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_obj", 2, 
    "_seq", seq.type,
    file_name_end))
  buffers[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_buffer", 
    "_obj", 1, 
    "_seq", seq.type,
    file_name_end))
  randoms[[i]] = readRDS(paste0(
    output_home, 
    "/scenario4_random", 
    file_name_end0))
  
  spacefills[[i]] = readRDS(paste0(
    output_home, 
    "/scenario4_spacefilling", 
    file_name_end0))
  seqmed.n1s[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_error", 1,
    "_obj", 1, 
    "_seq", seq.type,
    file_name_end0))
  seqmed.n2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_error", 2,
    "_obj", 1, 
    "_seq", seq.type,
    file_name_end0))
  seqmed.s1s[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_signal", 1,
    "_obj", 1, 
    "_seq", seq.type,
    file_name_end))
  seqmed.s2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario4_seqmed", 
    "_signal", 2,
    "_obj", 1, 
    "_seq", seq.type,
    file_name_end))
}

################################################################################
# make plots
################################################################################

# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
              error.var = nugget)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = nugget)
modelT = list(type = typeT, l = lT, signal.var = sigmasq, error.var = nugget)

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
                  aes(x = input, y = value, group = type, color = type, 
                      linetype = type)) + 
  geom_point() + 
  geom_path() +
  ylim(0, 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "P(HT|X, Y)", x = "Initial Data")
PPHT.plt
