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
sigmasq = 1
sigmasqs = c(1 - 1e-10, 1)
nuggets = c(1e-10, 1e-15)
nugget = NULL
nugget.sm = nugget # nugget for different signal variances and buffer
buffer = 0

# boxhill settings
nugget.bh = nugget # wouldn't run without
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

boxhills = readRDS(paste0(
  output_home, 
  "/scenario4_boxhill", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

qs = readRDS(paste0(
  output_home,
  "/scenario4_seqmed",
  "_obj", 2,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

buffers = readRDS(paste0(
  output_home,
  "/scenario4_buffer",
  "_obj", 1,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

randoms = readRDS(paste0(
  output_home, 
  "/scenario4_random", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

spacefills = readRDS(paste0(
  output_home, 
  "/scenario4_spacefilling", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

seqmeds.n1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s1 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s2 = readRDS(paste0(
  output_home,
  "/scenario4_seqmed", 
  "_obj", 1,
  "_signal", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

################################################################################
# make sequential EPPH plots
################################################################################
PPHs_seq = list()

# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
                 error.var = NULL)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = NULL)
modelT = list(type = typeT, l = lT, signal.var = sigmasq, error.var = NULL)

getPPHseq_scen3 = function(design, model0, model1, modelT){
  PPH0_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  PPH1_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  PPHT_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  # modelT$error.var = min(model0$error.var, model1$error.var, modelT$error.var)
  for(i in 1:length(as.vector(na.omit(design$y.new)))){
    y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:i])
    x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:i])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = rep(1 / 3, 3), 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1), 
        Evidence_gp(y.tmp, x.tmp, modelT)
      )
    )
    PPH0_seq[i] = PPHs.tmp[1]
    PPH1_seq[i] = PPHs.tmp[2]
    PPHT_seq[i] = PPHs.tmp[3]
  }
  if(length(PPH0_seq) < Nnew){
    PPH0_seq[(length(PPH0_seq) + 1):Nnew] = NA
    PPH1_seq[(length(PPH1_seq) + 1):Nnew] = NA
    PPHT_seq[(length(PPHT_seq) + 1):Nnew] = NA
  }
  return(data.frame(
    index = 1:Nnew, 
    PPH0 = PPH0_seq, 
    PPH1 = PPH1_seq, 
    PPHT = PPHT_seq
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
  PPH_seq.bh = getPPHseq_scen3(bh, model0, model1, modelT) # model0.bh, model1.bh, modelT)
  PPH_seq.q = getPPHseq_scen3(q, model0, model1, modelT) #  model0.other, model1.other, modelT)
  PPH_seq.b = getPPHseq_scen3(b, model0, model1, modelT) #  model0.other, model1.other, modelT)
  PPH_seq.r = getPPHseq_scen3(r, model0, model1, modelT) # model0.other, model1.other, modelT)
  PPH_seq.sf = getPPHseq_scen3(sf, model0, model1, modelT) # model0.other, model1.other, modelT)
  PPH_seq.n1 = getPPHseq_scen3(n1, model0, model1, modelT) # model0.n1, model1.n1, modelT)
  PPH_seq.n2 = getPPHseq_scen3(n2, model0, model1, modelT) # model0.n2, model1.n2, modelT)
  PPH_seq.s1 = getPPHseq_scen3(s1, model0, model1, modelT) # model0.s1, model1.s1, modelT)
  PPH_seq.s2 = getPPHseq_scen3(s2, model0, model1, modelT) # model0.s2, model1.s2, modelT)
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
PPHTmean_seq = aggregate(PPH_seq$PPHT, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPHTmean_seq) = c("index", "type", "value")
PPHTmean_seq$Hypothesis = "HT"

PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq, PPHTmean_seq)
ggplot(PPHmean_seq, aes(x = index, y = value, color = type, linetype = type)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  theme_bw() +
  ylim(0, 1)
