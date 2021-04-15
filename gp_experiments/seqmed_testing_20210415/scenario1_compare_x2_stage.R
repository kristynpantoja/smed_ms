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
# signalvar.type = 2 # 1 = phi0 sigmasq != 1, 2 = phi1 sigmasq != 1
input.type = 3 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
seq.type = 2 # 1 = fully sequential, 2 = stage-sequential 3x5

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
nugget.bh = 1e-10
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
# Scenario 1: Squared exponential vs. matern, true = matern
################################################################################
type01 = c("squaredexponential", "matern")
l01= c(0.01, 0.01) # SIM SETTING
# l01= c(0.1, 0.1) # DEMO SETTING

################################################################################
# models - BoxHill
model0.bh = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
                 error.var = nugget.bh)
model1.bh = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                 error.var = nugget.bh)

# models - q, buffer
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
# import matern functions
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario1_simulated_functions", 
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
  "/scenario1_boxhill", 
  "_input", input.type, 
  "_seed", rng.seed, 
  ".rds"
))

qs = readRDS(paste0(
  output_home,
  "/scenario1_seqmed",
  "_obj", 2,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

buffers = readRDS(paste0(
  output_home,
  "/scenario1_buffer",
  "_obj", 1,
  "_input", input.type,
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n1 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
  "_obj", 1,
  "_error", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.n2 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
  "_obj", 1,
  "_error", 2, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s1 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
  "_obj", 1,
  "_signal", 1, 
  "_input", input.type, 
  "_seq", seq.type,
  "_seed", rng.seed,
  ".rds"
))

seqmeds.s2 = readRDS(paste0(
  output_home,
  "/scenario1_seqmed", 
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
# 6 designs that aren't boxhill
idx = 1
designs = list(qs[[idx]], buffers[[idx]], 
               seqmeds.n1[[idx]], seqmeds.n2[[idx]], 
               seqmeds.s1[[idx]], seqmeds.s2[[idx]])
model0s = list(model0.other, model0.other, 
               model0.n1, model0.n2, model0.s1, model0.s2)
model1s = list(model1.other, model1.other,
               model1.n1, model1.n2, model1.s1, model1.s2)
design.names = c("q", "buffer", 
                 "nugget1", "nugget2", 
                 "signal1", "signal2")

# global settings
candidates = x_seq
batch.idx = 1
N2 = 2

# function for getting objective of x2 -- stages!
getx2 = function(design, objective.type, buffer, model0, model1){
  
  initD = design$x
  y = design$y
  
  initN = length(initD)
  if(length(y) != initN) stop("length of y does not match length of initial input data, initD")
  
  # check if any points in initD give Wasserstein distance of 0 (in which case we don't want to use it since 1/0 in q)
  # turns out, that's not necessary
  
  # old_initD = initD
  
  # posterior distribution of beta
  if(is.null(model0$error.var)){
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l))
  } else{
    Kinv0 = solve(getCov(initD, initD, model0$type, model0$l) + 
                    diag(rep(model0$error.var, initN)))
  }
  if(is.null(model1$error.var)){
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l))
  } else{
    Kinv1 = solve(getCov(initD, initD, model1$type, model1$l) + 
                    diag(rep(model1$error.var, initN)))
  }
  
  D = rep(NA, N2)
  D_ind = rep(NA, N2)
  if(batch.idx == 1){
    # -- Initialize 1st additional design point-- #
    w_candidates = sapply(candidates, FUN = function(x) WNgp(
      x, Kinv0, Kinv1, initD, y, model0, model1))
    w_opt = which.max(w_candidates)
    xopt = candidates[w_opt]
    is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))
  } else{
    is_x_max_in_initD = TRUE
  }
  if(is_x_max_in_initD){
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) obj_gp(
        x, NULL, 
        Kinv0, Kinv1, initD, y, p, k, alpha, buffer, objective.type, 
        model0, model1))
    if(all(f_min_candidates == Inf)){
      stop("SeqMEDgp_batch: all candidates result in objective function = Inf.")
    }
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D[1] = xnew
    D_ind[1] = f_opt
  } else{
    D[1] = xopt
    D_ind[1] = w_opt
  }
  
  if(N2 > 1){
    for(i in 2:N2){
      # Find f_opt: minimum of f_min
      f_min_candidates = sapply(
        candidates, 
        function(x) obj_gp(
          x, D[1:(i - 1)], 
          Kinv0, Kinv1, initD, y, p = 1, k = 4, alpha = 1, buffer, objective.type, 
          model0, model1))
      f_opt = which.min(f_min_candidates)
      xnew = candidates[f_opt]
      # Update set of design points (D) and plot new point
      D[i] = xnew
      D_ind[i] = f_opt
    }
  }
  
  # checks #
  if(xnew != design$x.new[2]) warning("x2 doesn't match!")
  if(f_opt != design$x.new.idx[2]) warning("x2.idx doesn't match!")
  
  return(list(
    objective = f_min_candidates, 
    x = xnew,
    x.idx = f_opt
  ))
}

# calculate x2 objective
objectives = matrix(NA, nrow = length(x_seq), ncol = length(designs))
x_vec = rep(NA, length(designs))
x_idx_vec = rep(NA, length(designs))
for(i in 1:length(designs)){
  # objective.type
  if(i == 1){
    objective.type = 2
  } else{
    objective.type = 1
  }
  # buffer
  if(i == 2){
    buffer = 1e-15
  } else{
    buffer = 0
  }
  res = getx2(designs[[i]], objective.type, buffer, model0s[[i]], model1s[[i]])
  objectives[, i] = res$objective
  x_vec[i] = res$x
  x_idx_vec[i] = res$x.idx
}

design.names = c("q", "buffer", "nugget1", "nugget2", 
                 "signal1", "signal2")
design.levels = c("nugget1", "nugget2", "q", 
                 "signal1", "signal2", "buffer")

obj.data = data.frame(melt(objectives))
names(obj.data) = c("index", "type", "logObjective")
obj.data$logObjective = log(obj.data$logObjective)
obj.data$type = factor(obj.data$type, levels = 1:length(designs), 
                       labels = design.names) # map design numbers to names
# reorder
obj.data$type = factor(obj.data$type, levels = design.levels)
obj.data$x = rep(x_seq, length(designs))

ggplot(obj.data, aes(x = x, y = logObjective, color = )) + 
  facet_wrap(~type) + 
  geom_path()
