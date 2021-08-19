################################################################################
# last updated: 08/17/2021
# purpose: to test seqmedvs for dimensions 2,3

dimT = 2 # 2, 3

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lmvs/outputs"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
# source(paste(functions_home, "/MMED.R", sep = ""))
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
numSeq = 9
Nin = 5
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = -1
xmax = 1
dimX = 3
numCandidates = 5000
x_seq = seq(from = xmin, to = xmax, length.out = floor((numCandidates)^(1 / 3)))
candidates = expand.grid(x_seq, x_seq, x_seq)
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
  fT = function(x) betaT %*% x[indicesT]
} else if(dimT == 3){
  betaT = mu_full[indices1]
  indicesT = indices0
  fT = function(x) betaT %*% x[indicesT]
}

# generate seqmeds
registerDoRNG(rng.seed)
seqmed_list = foreach(i = 1:numSims) %dorng% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  SeqMED()
  
}
saveRDS(seqmed_list, paste(output_home, "/seqmed/dim2_seqmed_simulations", 
                           "_numSeq", numSeq, 
                           "_seqN", seqN,
                           "_numSims", numSims,
                           ".rds", sep = ""))



