################################################################################
# last updated: 08/18/21
# purpose: to create a list of seqmed simulations
# dimT = 2:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2)
# dimT = 3:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2, 3)

dimT = 2 # 2, 3

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lmvs/experimenting"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/construct_design_matrix.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
# source(paste(functions_dir, "/MMED.R", sep = ""))
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
# source(paste(functions_dir, "/boxhill_vs.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

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
Nin = 5
numSeq = 27
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = -1
xmax = 1
dimX = 3
numCandidates = 5000
x_seq = seq(from = xmin, to = xmax, length.out = floor((numCandidates)^(1 / 3)))
candidates = as.matrix(expand.grid(x_seq, x_seq, x_seq))
sigmasq = 0.3

# hypothesis settings
type01 = c(2, 3)
mu_full = c(0.5, 0.5, 0.5) #
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
sigmasq01 = 0.5
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
model0 = list(
  indices = indices0, beta.mean = mu0, beta.var = V0)
model1 = list(
  indices = indices1, beta.mean = mu1, beta.var = V1)

# seqmed settings
p = 3
k = 4 * p

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# scenario settings
################################################################################
if(dimT == 2){
  betaT = mu_full[indices0]
  indicesT = indices0
  fT = function(x) x[, indicesT, drop = FALSE] %*% betaT
} else if(dimT == 3){
  betaT = mu_full[indices1]
  indicesT = indices0
  fT = function(x) x[, indicesT, drop = FALSE] %*% betaT
}

# generate seqmeds -- new
registerDoRNG(rng.seed)
initD = matrix(runif(n = dimX * Nin, min = xmin, max = xmax), nrow = Nin, ncol = dimX)
inity = as.vector(simulateYvs(initD[ , indicesT], Nin, betaT, sigmasq, 1, seed = seed))
sm.res = SeqMEDvs(
  y.in = inity, x.in = initD, model0 = model0, model1 = model1, 
  error.var = sigmasq, candidates = candidates, true.function = fT, 
  true.indices = indicesT, dimX = dimX, k = k, xmin = xmin, xmax = xmax, p = p, 
  numSeq = numSeq, seqN = seqN, alpha_seq = 1, prints = FALSE, seed = NULL
  # mean_beta_full = mu_full, beta_true = betaT, indices_true = indicesT, 
  # indices0 = indices0, indices1 = indices1, mean_beta0 = mu0, mean_beta1 = mu1,
  # var_e = sigmasq, var_beta = sigmasq01, var_beta0 = V0, var_beta1 = V1, 
  # xmin = xmin, xmax = xmax, numCandidates = numCandidates, k = k, p = p,
  # initD = initD, inity = inity, numSeq = numSeq, N_seq = seqN, 
  # alpha_seq = NULL, buffer_seq = 0, candidates = candidates, wasserstein0 = 1, 
  # genCandidates = 1, seed = 1
)



# generate seqmeds -- old
# registerDoRNG(rng.seed)
# f0 = function(x) mu0 %*% x[indices0]
# f1 = function(x) mu1 %*% x[indices1]
# # initial design
# initD = matrix(runif(n = dimX * Nin, min = xmin, max = xmax), nrow = Nin, ncol = dimX)
# inity = as.vector(simulateYvs(initD[ , indicesT], Nin, betaT, sigmasq, 1, seed = seed))
# sm.res = generate_SMMEDvs(
#   mean_beta_full = mu_full, beta_true = betaT, indices_true = indicesT, 
#   indices0 = indices0, indices1 = indices1, mean_beta0 = mu0, mean_beta1 = mu1,
#   var_e = sigmasq, var_beta = sigmasq01, var_beta0 = V0, var_beta1 = V1, 
#   xmin = xmin, xmax = xmax, numCandidates = numCandidates, k = k, p = p,
#   initD = initD, inity = inity, numSeq = numSeq, N_seq = seqN, 
#   alpha_seq = NULL, buffer_seq = 0, candidates = candidates, wasserstein0 = 1, 
#   genCandidates = 1, seed = 1)
# saveRDS(sm.res, paste(output_dir, "/seqmed.rds", sep = ""))
sm.res.old = readRDS(paste(output_dir, "/seqmed.rds", sep = ""))
all.equal(rbind(sm.res$x.in, sm.res$x.new), sm.res.old$D)


# mean_beta_full = mu_full
# beta_true = betaT
# indices_true = indicesT
# # indices0 = indices0, indices1 = indices1, 
# mean_beta0 = mu0
# mean_beta1 = mu1
# var_e = sigmasq
# var_beta = sigmasq01
# var_beta0 = V0
# var_beta1 = V1
# # xmin = xmin, xmax = xmax, numCandidates = numCandidates, k = k, p = p,
# # initD = initD, inity = inity, numSeq = numSeq, 
# N_seq = seqN
# alpha_seq = NULL
# buffer_seq = 0
# candidates = candidates
# wasserstein0 = 1
# genCandidates = 1
# seed = 1

hist(sm.res$x.new[, 3])




