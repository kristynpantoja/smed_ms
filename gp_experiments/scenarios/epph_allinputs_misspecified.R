################################################################################
# last updated: 04/28/2021
# purpose: to test seqmedgp for scenario 3:
#   squared exponential vs. another squared exponential,
#   where the true function is matern

scenario = 3 # scenarios: 3, 4, 5, 6
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5

################################################################################
# Sources/Libraries
################################################################################
output_home = paste0("gp_experiments/scenarios/scenario", scenario, "/outputs")
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
sigmasq_err = 1e-10

# SeqMED settings
sigmasqs = c(1 - 1e-10, 1)
if(scenario %in% c(3, 4)){
  nuggets = c(1e-10, 1e-15)
} else if(scenario %in% c(5, 6)){
  nuggets = c(1e-5, 1e-10)
}

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
# Scenario settings
################################################################################
if(scenario == 3){
  type01 = c("squaredexponential", "squaredexponential")
  typeT = "matern"
  l01= c(0.005, 0.01)
  lT = 0.01
} else if(scenario == 4){
  type01 = c("matern", "squaredexponential")
  typeT = "periodic"
  l01= c(0.01, 0.01)
  lT = 0.01
} else if(scenario == 5){
  type01 = c("matern", "periodic")
  typeT = "squaredexponential"
  l01= c(0.01, 0.01)
  lT = 0.01
} else if(scenario == 6){
  type01 = c("squaredexponential", "periodic")
  typeT = "matern"
  l01= c(0.01, 0.01)
  lT = 0.01
}

################################################################################
# import matern functions
filename_append = ""
if(!is.null(sigmasq_err)){
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_err), "-")[[1]][2])
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
  seqmed.s1s[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
    "_signal", 1,
    "_seq", seq.type,
    filename_append.tmp))
  seqmed.s2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
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
    "/scenario", scenario, "_random", 
    filename_append.tmp))
  spacefills[[i]] = readRDS(paste0(
    output_home, 
    "/scenario", scenario, "_spacefilling", 
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
    "/scenario", scenario, "_seqmed", 
    "_error", 1,
    "_seq", seq.type,
    filename_append.tmp))
  seqmed.n2s[[i]] = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_seqmed", 
    "_error", 2,
    "_seq", seq.type,
    filename_append.tmp))
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
  index = length(as.vector(na.omit(design$y.new)))
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
  return(data.frame(
    index = index, 
    "H0" = PPHs.tmp[1], "H1" = PPHs.tmp[2], "HT" = PPHs.tmp[3]))
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
    PPH.bh = getPPH(bh, model0, model1, modelT)
    PPH.q = getPPH(q, model0, model1, modelT)
    PPH.b = getPPH(b, model0, model1, modelT)
    PPH.r = getPPH(r, model0, model1, modelT)
    PPH.sf = getPPH(sf, model0, model1, modelT)
    PPH.n1 = getPPH(n1, model0, model1, modelT)
    PPH.n2 = getPPH(n2, model0, model1, modelT)
    PPH.s1 = getPPH(s1, model0, model1, modelT)
    PPH.s2 = getPPH(s2, model0, model1, modelT)
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

### average over all, regardless of index of last point

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
  labs(y = "E[P(HT|X, Y)|X]", x = "Initial Data")
PPHT.plt


# ### average over largest index of last point
# 
# # data.max.index = dplyr::filter(PPH, index == max(PPH$index))
# PPH0mean2 = aggregate(PPH$H0, 
#                      by = list(PPH$type, PPH$input, PPH$index), 
#                      FUN = function(x) mean(x, na.rm = TRUE))
# names(PPH0mean2) = c("type", "input", "index", "value")
# PPH0mean2$Hypothesis = "H0"
# PPH1mean2 = aggregate(PPH$H1, 
#                       by = list(PPH$type, PPH$input, PPH$index), 
#                      FUN = function(x) mean(x, na.rm = TRUE))
# names(PPH1mean2) = c("type", "input", "index", "value")
# PPH1mean2$Hypothesis = "H1"
# PPHTmean2 = aggregate(PPH$HT, 
#                       by = list(PPH$type, PPH$input, PPH$index),  
#                      FUN = function(x) mean(x, na.rm = TRUE))
# names(PPHTmean2) = c("type", "input", "index", "value")
# PPHTmean2$Hypothesis = "HT"
# 
# PPHmean2 = rbind(PPH0mean2, PPH1mean2, PPHTmean2)
# PPHmean2$type = factor(PPHmean2$type)
# PPHmean2$Hypothesis = factor(PPHmean2$Hypothesis)
# PPHmean2$input = factor(PPHmean2$input, 
#                        labels = c("extrapolation", "inc spread", "even coverage"))
# 
# PPHT.plt2 = ggplot(dplyr::filter(PPHmean2, Hypothesis == "HT"), 
#                   aes(x = input, y = value, group = type, color = type, 
#                       linetype = type)) + 
#   geom_point() + 
#   geom_line() +
#   ylim(0, 1) + 
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   labs(y = "E[P(HT|X, Y)|X]", x = "Initial Data")
# PPHT.plt2
