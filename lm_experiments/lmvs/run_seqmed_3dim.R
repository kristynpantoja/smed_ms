
# dimT = 2:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2)
# dimT = 3:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2, 3)

dimT = 3 # 2, 3
sigmasq = 0.1 # 0.1, 0.3

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lmvs/outputs"
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
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

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

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
Nin = 1 #1, 5
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

# hypothesis settings
mu_full = rep(0.5, 3)
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
sigmasq01 = 0.25 # 0.5
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
model0 = list(
  indices = indices0, beta.mean = mu0, beta.var = V0)
model1 = list(
  indices = indices1, beta.mean = mu1, beta.var = V1)

# seqmed settings
p = dimX
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
  indicesT = indices1
  fT = function(x) x[, indicesT, drop = FALSE] %*% betaT
}

################################################################################
# run simulations
################################################################################

# generate seqmeds
registerDoRNG(rng.seed)
seqmed_list = foreach(i = 1:numSims) %dorng% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  initD = matrix(runif(
    n = dimX * Nin, min = xmin, max = xmax), nrow = Nin, ncol = dimX)
  inity = simulateY_frommultivarfunction(
    x = initD[, indicesT, drop = FALSE], true.function = fT, 
    error.var = sigmasq)
  SeqMEDvs(
    y.in = inity, x.in = initD, model0 = model0, model1 = model1, 
    error.var = sigmasq, candidates = candidates, true.function = fT, 
    true.indices = indicesT, dimX = dimX, k = k, xmin = xmin, xmax = xmax, 
    p = p, numSeq = numSeq, seqN = seqN, alpha_seq = 1, prints = FALSE)
}
saveRDS(seqmed_list, paste0(
  output_dir, "/3dim", 
  "_dim", dimT, 
  "_seqmed", 
  "_Nttl", Nttl,
  "_Nin", Nin,
  "_numSeq", numSeq,
  "_seqN", seqN,
  "_seed", rng.seed,
  "_noise", strsplit(as.character(sigmasq), split = "\\.")[[1]][2],
  ".rds"))


