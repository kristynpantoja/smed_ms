################################################################################
# last updated: 04/07/2021
# purpose: to test seqmedgp for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern
# trying out some (not necessarily MED) designs

################################################################################
# Sources/Libraries
################################################################################
output_home = "run_designs/gp_experiments"
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
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)
library(doRNG)

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
numSims = 100
seed = 12
Nin = 6
numSeq = 15
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)

# SeqMED settings
signalSM = 1
nuggetSM = NULL
buffer = 1e-20 # NONZERO #######################################################

# boxhill settings
prior_probs = rep(1 / 2, 2)
nuggetBH = 1e-10

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

# generate matern functions ####################################################
set.seed(seed)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values

# bh settings
model0 = list(type = type01[1], l = l01[1], signal.var = 1, error.var = nuggetSM)
model1 = list(type = type01[2], l = l01[2], signal.var = 1, error.var = nuggetSM)

################################################################################
# generate seqmeds #############################################################
################################################################################

################################################################################
# try different objective.type - SIM SETTINGS ##################################
################################################################################

# input set
x_input = x_in1
x_input_idx = x_in1_idx
seqmed_list = list()
# index
i = 1
y_seq = y_seq_mat[ , i]
y_input = y_seq[x_input_idx]

# calculate posterior probs given initial data
getHypothesesPosteriors(
  prior.probs = prior_probs, 
  evidences = c(
    Evidence_gp(y_input, x_input, model0, nuggetBH),
    Evidence_gp(y_input, x_input, model1, nuggetBH)
  )
)

# seqmeds
types = c("BH", "MED", "q", "BatchMED", "Batch q")
obj.seq = c(1, 2)
seqmed_list = list()
for(i in 1:length(obj.seq)){
  numSeq = 15 # FULL SEQUENTIAL ################################################
  seqN = 1 # FULL SEQUENTIAL ###################################################
  Nnew = numSeq * seqN # FULL SEQUENTIAL #######################################
  seqmed_list[[i]] = SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, model0 = model0, model1 = model1, 
    numSeq = numSeq, seqN = seqN, prints = TRUE, buffer = buffer, 
    objective.type = obj.seq[i] # FOR OBJECTIVE DEMO ###########################
    , seed = 1234
  )
}
for(i in 1:length(obj.seq)){
  idx = i + 2
  numSeq = 3 # BATCH SEQUENTIAL ################################################
  seqN = 5 # BATCH SEQUENTIAL ##################################################
  Nnew = numSeq * seqN # BATCH SEQUENTIAL ######################################
  seqmed_list[[idx]] = SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, model0 = model0, model1 = model1,
    numSeq = numSeq, seqN = seqN, prints = TRUE, buffer = buffer, 
    objective.type = obj.seq[i] # FOR OBJECTIVE DEMO ###########################
    , seed = 1234
  )
}
# boxhill
bh = BHgp_m2(
  y_input, x_input, x_input_idx, prior_probs, model0, model1, Nnew, 
  x_seq, y_seq, nuggetBH)
# bh$x.new

# par(mfrow = c(4, 1))
x.new.mat = matrix(NA, nrow = Nnew, ncol = length(seqmed_list))
for(i in 1:length(seqmed_list)){
  x.in.tmp = seqmed_list[[i]]$x
  x.new.tmp = seqmed_list[[i]]$x.new
  x.new.mat[, i] = x.new.tmp
  # plot(x.in.tmp, y = rep(0, length(x.in.tmp)), 
  #      ylim = c(-0.01, 0.02), xlim = c(xmin, xmax),
  #      xlab = "", ylab = "")
  # points(x.new.tmp, y = rep(0.01, length(x.new.tmp)), col = 2)
  # print(x.new.tmp)
}
x.new.mat = cbind(bh$x.new, x.new.mat)

# par(mfrow = c(1, 1))
data.gg = data.frame(
  Index = as.character(rep(1:Nnew, length(seqmed_list) + 1)), 
  Sim = rep(1:(length(seqmed_list) + 1), each = Nnew), 
  Type = factor(rep(types, each = Nnew), levels = types), 
  SeqMEDgp = as.vector(x.new.mat)
)
data.gg0 = data.frame(
  Sim = rep(1:(length(seqmed_list) + 1), each = Nin), 
  Type = factor(rep(types, each = Nin), levels = types), 
  Input = rep(x_input, (length(seqmed_list) + 1))
)
ggplot() + 
  geom_point(data = data.gg0, 
             mapping = aes(x = Input, y = Type)) +
  geom_point(data = data.gg, 
             mapping = aes(x = SeqMEDgp, y = Type, color = Type), 
             inherit.aes = FALSE) + 
  geom_text(data = data.gg, aes(x = SeqMEDgp, y = Type, label = Index), 
            vjust = -0.5 * as.numeric(paste(data.gg$Index)), size = 1.75) +
  xlim(c(xmin, xmax))

# plot the function
# ggplot(data.frame(x = x_seq, y = y_seq), aes(x = x, y = y)) + 
#   geom_path()

################################################################################
# posterior probabilities ######################################################
################################################################################


PPHs = list()
PPHs_seq = list()

# PPHs for BH
y.tmp = c(bh$y, as.vector(na.omit(bh$y.new)))
x.tmp = c(bh$x, as.vector(na.omit(bh$x.new)))
PPHs[[1]] = getHypothesesPosteriors(
  prior.probs = prior_probs, 
  evidences = c(
    Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
    Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
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
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
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

################################################################################
# objectives for 1st point #####################################################
################################################################################
# this is the same for MED, q, BatchMED, Batch q
initD = x_input
y = y_input
candidates = x_seq

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

# -- Initialize 1st additional design point-- #
w_candidates = sapply(candidates, FUN = function(x) WNgp(
  x, Kinv0, Kinv1, initD, y, model0, model1))
w_opt = which.max(w_candidates)
xopt = candidates[w_opt]
is_x_max_in_initD = any(sapply(initD, function(x) x == xopt))

ggplot(data.frame(x = candidates, y = w_candidates), aes(x = x, y = y)) + 
  geom_path() + 
  geom_point(data = data.frame(x = x_input, y = 0), 
             mapping = aes(x = x, y = y), color = "red",
             inherit.aes = FALSE) + 
  geom_point(data = data.frame(x = xopt, y = w_candidates[w_opt]), 
             mapping = aes(x = x, y = y), color = "green",
             inherit.aes = FALSE)

################################################################################
# objectives for 2nd point #####################################################
################################################################################
# where the methods differ :0

objectives = matrix(NA, nrow = length(x_seq), ncol = length(seqmed_list))
xopt_idx_objectives = rep(NA, length(seqmed_list))
xopt_objectives = rep(NA, length(seqmed_list))

for(i in 1:length(obj.seq)){
  seqmed.tmp = seqmed_list[[i]]
  objective.type = obj.seq[i]
  batch.idx = 2
  
  initD = c(seqmed.tmp$x, seqmed.tmp$x.new[1])
  y = c(seqmed.tmp$y, seqmed.tmp$y.new[1])
  
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
  
    # Find f_opt: minimum of f_min
    f_min_candidates = sapply(
      candidates, 
      function(x) obj_gp(
        x, NULL, 
        Kinv0, Kinv1, initD, y, p = 1, k = 4, alpha = 1, buffer, objective.type, 
        model0, model1))
    if(all(f_min_candidates == Inf)){
      stop("SeqMEDgp_batch: all candidates result in objective function = Inf.")
    }
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
  
  # assign #
  objectives[, i] = f_min_candidates
  xopt_idx_objectives[i] = xnew
  xopt_objectives[i] = f_opt
  
  # checks #
  print(xnew == seqmed.tmp$x.new[2])
  print(f_opt == seqmed.tmp$x.new.idx[2])
}
for(i in 1:length(obj.seq)){
  idx = i + 2
  seqmed.tmp = seqmed_list[[idx]]
  objective.type = obj.seq[i]
  batch.idx = 1
  N2 = 2
  
  initD = seqmed.tmp$x
  y = seqmed.tmp$y
  
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
  
  # assign #
  objectives[, idx] = f_min_candidates
  xopt_idx_objectives[idx] = xnew
  xopt_objectives[idx] = f_opt
  
  # checks #
  print(xnew == seqmed.tmp$x.new[2])
  print(f_opt == seqmed.tmp$x.new.idx[2])
}

colnames(objectives) = types[-1]
objectives.df = melt(objectives)
objectives.df = as.data.frame(objectives.df)
names(objectives.df) = c("Index", "Type", "value")
objectives.df$x = rep(x_seq, ncol(objectives))
objectives.df$value = log(objectives.df$value)

ggplot(objectives.df, aes(x = x, y = value, color = Type)) + 
  facet_wrap(~Type, scales = "free_y") + geom_path()

