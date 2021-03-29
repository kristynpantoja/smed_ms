################################################################################
# last updated: 03/28/2021
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
nuggetSM = 1e-10

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
model0 = list(type = type01[1], l = l01[1])
model1 = list(type = type01[2], l = l01[2])

################################################################################
# generate seqmeds #############################################################
################################################################################

################################################################################
# try different buffers ########################################################
################################################################################

# input set
x_input = x_in3
x_input_idx = x_in3_idx
seqmed_list = list()
# index
i = 1
y_seq = y_seq_mat[ , i]
y_input = y_seq[x_input_idx]
# seqmed
buffer.seq = c(1e-20, 1e-15, 1e-10)
seqmed_list = list()
for(i in 1:length(buffer.seq)){
  seqmed_list[[i]] = SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, nugget = nuggetSM, type = type01, l = l01,
    numSeq = numSeq, seqN = seqN, prints = TRUE, buffer = buffer.seq[i]
    , seed = 1234 # DELETE THIS ARGUMENT LATER.
  )
}
# y0 = y_input
# x0 = x_input
# x0.idx = x_input_idx
# candidates = x_seq
# function.values = y_seq
# nugget = nuggetSM
# type = type01
# l = l01
# error.var = 1
# xmin = 0
# xmax = 1
# k = 4
# p = 1
# numSeq
# seqN
# alpha.seq = 1
# buffer = buffer.seq[i]
# objective.type = 1
# init.as.stage = FALSE
# prints = TRUE
# seed = NULL


par(mfrow = c(5, 1))
x.new.mat = matrix(NA, nrow = Nnew, ncol = length(buffer.seq))
for(i in 1:length(buffer.seq)){
  x.in.tmp = seqmed_list[[i]]$x
  x.new.tmp = seqmed_list[[i]]$x.new
  x.new.mat[, i] = x.new.tmp
  plot(x.in.tmp, y = rep(0, length(x.in.tmp)), 
       ylim = c(-0.01, 0.02), xlim = c(xmin, xmax),
       xlab = "", ylab = "")
  points(x.new.tmp, y = rep(0.01, length(x.new.tmp)), col = 2)
  print(x.new.tmp)
}

par(mfrow = c(1, 1))
data.gg = data.frame(
  Sim = rep(1:length(buffer.seq), each = Nnew), 
  Buffer = factor(rep(buffer.seq, each = Nnew), levels = buffer.seq), 
  SeqMEDgp = as.vector(x.new.mat)
  )
data.gg0 = data.frame(
  Sim = rep(1:length(buffer.seq), each = Nin), 
  Buffer = factor(rep(buffer.seq, each = Nin), levels = buffer.seq), 
  Input = rep(x_input, length(buffer.seq))
)
ggplot() + 
  geom_point(data = data.gg0, 
             mapping = aes(x = Input, y = Buffer)) +
  geom_point(data = data.gg, 
             mapping = aes(x = SeqMEDgp, y = Buffer, color = Buffer), 
             inherit.aes = FALSE) + 
  xlim(c(xmin, xmax))

# plot the function
ggplot(data.frame(x = x_seq, y = y_seq), aes(x = x, y = y)) + 
  geom_path()

################################################################################
# try different objective.type #################################################
################################################################################

# input set
x_input = x_in1
x_input_idx = x_in1_idx
seqmed_list = list()
# index
i = 1
y_seq = y_seq_mat[ , i]
y_input = y_seq[x_input_idx]
# seqmed
buffer = 1e-20
types = c("MED", "q", "BatchMED", "Batch q")
obj.seq = c(1, 2)
seqmed_list = list()
for(i in 1:length(obj.seq)){
  seqmed_list[[i]] = SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, nugget = nuggetSM, type = type01, l = l01,
    numSeq = numSeq, seqN = seqN, prints = TRUE, buffer = buffer, 
    objective.type = obj.seq[i]
    , seed = 1234 # DELETE THIS ARGUMENT LATER.
  )
}
for(i in 1:length(obj.seq)){
  idx = i + 2
  seqmed_list[[idx]] = SeqMEDgp(
    y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
    function.values = y_seq, nugget = nuggetSM, type = type01, l = l01,
    numSeq = 3, seqN = 5, prints = TRUE, buffer = buffer, 
    objective.type = obj.seq[i]
    , seed = 1234 # DELETE THIS ARGUMENT LATER.
  )
}

par(mfrow = c(4, 1))
x.new.mat = matrix(NA, nrow = Nnew, ncol = length(seqmed_list))
for(i in 1:length(seqmed_list)){
  x.in.tmp = seqmed_list[[i]]$x
  x.new.tmp = seqmed_list[[i]]$x.new
  x.new.mat[, i] = x.new.tmp
  plot(x.in.tmp, y = rep(0, length(x.in.tmp)), 
       ylim = c(-0.01, 0.02), xlim = c(xmin, xmax),
       xlab = "", ylab = "")
  points(x.new.tmp, y = rep(0.01, length(x.new.tmp)), col = 2)
  print(x.new.tmp)
}

par(mfrow = c(1, 1))
data.gg = data.frame(
  Sim = rep(1:length(seqmed_list), each = Nnew), 
  Type = factor(rep(types, each = Nnew), levels = types), 
  SeqMEDgp = as.vector(x.new.mat)
)
data.gg0 = data.frame(
  Sim = rep(1:length(seqmed_list), each = Nin), 
  Type = factor(rep(types, each = Nin), levels = types), 
  Input = rep(x_input, length(seqmed_list))
)
ggplot() + 
  geom_point(data = data.gg0, 
             mapping = aes(x = Input, y = Type)) +
  geom_point(data = data.gg, 
             mapping = aes(x = SeqMEDgp, y = Type, color = Type), 
             inherit.aes = FALSE) + 
  xlim(c(xmin, xmax))

# plot the function
ggplot(data.frame(x = x_seq, y = y_seq), aes(x = x, y = y)) + 
  geom_path()
