################################################################################
# last updated: 05/19/2021
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
sigmasq_err = 1e-10

# SeqMED settings
sigmasqs = c(1 - 1e-10, 1)
nuggets = c(1e-5, 1e-10)

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
# Scenario 2: Matern vs. periodic, true = periodic
################################################################################
type01 = c("matern", "periodic")
typeT = type01[2]
l01= c(0.01, 0.01)
lT = l01[2]

################################################################################
# import periodic functions
filename_append = ""
if(!is.null(sigmasq_err)){
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_err), "-")[[1]][2])
}
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario2_simulated_functions", filename_append,
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
  # filename_append.tmp for boxhills, buffers, qs, signal seqmeds
  filename_append.tmp = filename_append
  if(!is.null(nugget)){
    filename_append.tmp = paste0(
      filename_append.tmp, 
      "_nugget", strsplit(as.character(nugget), "-")[[1]][2])
  }
  filename_append.tmp = paste0(
    filename_append.tmp, 
    "_input", i, 
    "_seed", rng.seed,
    ".rds")
  boxhills[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_boxhill", 
    filename_append.tmp))
  buffers[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_buffer", 
    "_seq", seq.type,
    filename_append.tmp))
  qs[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_q",
    "_seq", seq.type,
    filename_append.tmp))
  seqmed.s1s[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_signal", 1,
    "_seq", seq.type,
    filename_append.tmp))
  seqmed.s2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_signal", 2,
    "_seq", seq.type,
    filename_append.tmp))
  
  # filename_append.tmp for random, space-filling
  filename_append.tmp = filename_append
  filename_append.tmp = paste0(
    filename_append.tmp, 
    "_input", i, 
    "_seed", rng.seed,
    ".rds")
  randoms[[i]] = readRDS(paste0(
    output_home, 
    "/scenario2_random", 
    filename_append.tmp))
  spacefills[[i]] = readRDS(paste0(
    output_home, 
    "/scenario2_spacefilling", 
    filename_append.tmp))
  
  # filename_append.tmp for error seqmeds
  filename_append.tmp = filename_append
  nuggets_vals = strsplit(as.character(nuggets), "-")
  nuggets_vals = paste0(nuggets_vals[[1]][2], nuggets_vals[[2]][2])
  filename_append.tmp = paste0(
    filename_append.tmp,
    "_nuggets", nuggets_vals)
  filename_append.tmp = paste0(
    filename_append.tmp, 
    "_input", i, 
    "_seed", rng.seed,
    ".rds"
  )
  seqmed.n1s[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_error", 1,
    "_seq", seq.type,
    filename_append.tmp))
  seqmed.n2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario2_seqmed", 
    "_error", 2,
    "_seq", seq.type,
    filename_append.tmp))
}

################################################################################
# make plots
################################################################################
PPHs_seq = list()

# input set
bh.in = boxhills[[input.type]]
q.in = qs[[input.type]]
buf.in = buffers[[input.type]]
ran.in = randoms[[input.type]]
sf.in = spacefills[[input.type]]
n1.in = seqmed.n1s[[input.type]]
n2.in = seqmed.n2s[[input.type]]
s1.in = seqmed.s1s[[input.type]]
s2.in = seqmed.s2s[[input.type]]

# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
              error.var = nugget)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              error.var = nugget)

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
  bh = bh.in[[j]]
  q = q.in[[j]]
  b = buf.in[[j]]
  r = ran.in[[j]]
  sf = sf.in[[j]]
  n1 = n1.in[[j]]
  n2 = n2.in[[j]]
  s1 = s1.in[[j]]
  s2 = s2.in[[j]]
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
