
# scenario 1:
#   linear vs. quadratic,
#   where the true function is quadratic
# scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic
rm(list = ls())

scenario = 1 # 1, 2

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lm/outputs"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMED.R", sep = ""))
source(paste(functions_dir, "/SeqMED_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/construct_design_matrix.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
source(paste(functions_dir, "/MMED.R", sep = ""))
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill deisign
source(paste(functions_dir, "/boxhill.R", sep = ""))

# set up parallelization
# library(foreach)
# library(future)
# library(doFuture)
# library(parallel)
# registerDoFuture()
# nworkers = detectCores() / 2
# plan(multisession, workers = nworkers)

# library(rngtools)
# library(doRNG)
rng.seed = 123 # 123
# registerDoRNG(rng.seed)
set.seed(123)

library(microbenchmark)

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 500 #
numSeq = 100 # 12, 100, 15 for timing
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
if(scenario == 1){
  if(numSeq == 100){
    sigmasq = 0.28
  } else if(numSeq == 12 | numSeq == 15){
    sigmasq = 0.04
  }
} else if(scenario == 2){
  if(numSeq == 100){
    sigmasq = 0.21
  } else if(numSeq == 1 | numSeq == 152){
    sigmasq = 0.038
  }
}

# shared settings
type01 = c(2, 3)
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
desX0 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x))
}
desX1 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x, x^2))
}
model0 = list(designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(designMat = desX1, beta.mean = mu1, beta.var = V1)

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenarios
################################################################################
if(scenario == 1){
  betaT = c(-0.2, -0.4, 0.4)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
} else if(scenario == 2){
  betaT = c(0, -0.75, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
}

################################################################################
# run simulations
################################################################################

# generate seqmeds
i = 1

# start.time = Sys.time()
# seqmed = SeqMED(
#   y.in = NULL, x.in = NULL, true.function = fT,
#   model0 = model0, model1 = model1, 
#   error.var = sigmasq, xmin = xmin, xmax = xmax,
#   candidates = candidates, numSeq = numSeq, seqN = seqN, 
#   alpha_seq = 1)
# end.time = Sys.time()
# time.diff = difftime(end.time, start.time, units = "secs")
# time.diff

microbenchmark(
  seqmed = SeqMED(
    y.in = NULL, x.in = NULL, true.function = fT,
    model0 = model0, model1 = model1, 
    error.var = sigmasq, xmin = xmin, xmax = xmax,
    candidates = candidates, numSeq = numSeq, seqN = seqN, 
    alpha_seq = 1), 
  times = 10
)




