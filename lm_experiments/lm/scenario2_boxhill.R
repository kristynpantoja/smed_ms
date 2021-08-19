################################################################################
# last updated: 12/09/20
# purpose: to create a list of boxhill simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic

################################################################################
# Sources/Libraries
################################################################################
output_home = "lm_experiments/lm/outputs"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMED.R", sep = ""))
source(paste(functions_home, "/SeqMED_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

# for generating initial data
source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))

# for box-hill deisign
source(paste(functions_home, "/boxhill.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(knitr)

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(doRNG)
registerDoRNG(1995)

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100
numSeq = 100
seqN = 1
N = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
sigmasq = 0.1

# shared settings
type01 = c(2, 3)
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
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
# Scenario 2: True function is cubic
################################################################################
betaT = c(0, -0.75, 0, 1)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3

################################################################################
# run simulations
################################################################################

# generate boxhills
bh_list = foreach(i = 1:numSims) %dorng% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  if(MMEDinputdata){
    N.new = (numSeq - 1) * seqN
    seqmed.res = SeqMED(
      D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
      beta.mean0 = mu0, beta.mean1 = mu1, beta.var0 = V0, beta.var1 = V1, 
      error.var = sigmasq, f0 = f0, f1 = f1, type = type01, xmin = xmin, xmax = xmax, 
      candidates = candidates, numSeq = 1, seqN = seqN
    )
    x_input = seqmed.res$D
    y_input = seqmed.res$y
    bh.res = BH_m2(y_input, x_input, prior_probs, model0, model1, N.new, 
                   candidates, fT, sigmasq)
  } else{
    bh.res = BH_m2(NULL, NULL, prior_probs, model0, model1, N, 
                   candidates, fT, sigmasq)
  }
  bh.res
}
saveRDS(bh_list, 
        paste(output_home, "/boxhill/scenario2_boxhill_simulations", 
              "_N", N, 
              "_MMEDinput", as.numeric(MMEDinputdata), 
              "_numSims", numSims, 
              ".rds", sep = ""))
