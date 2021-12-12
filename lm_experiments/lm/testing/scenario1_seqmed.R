################################################################################
# last updated: 10/09/2021
# purpose: to create a list of seqmed simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is quadratic
rm(list = ls())

scenario = 1
sigmasq = 0.2 # 0.1, 0.05, 0.025
numSeq = 12 #100, 36, 12
numSims = 100 #100
alpha = 0
sequential_alpha = FALSE
hybrid_alpha = FALSE
if(sequential_alpha) alpha_seq = seq(0, alpha, numSeq)
if(hybrid_alpha) alpha_seq = c(rep(0, numSeq / 2), rep(alpha, numSeq / 2))

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lm/testing"
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

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))

# for D-optimal design
library(AlgDesign)

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggpubr)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

library(doRNG)
registerDoRNG(1995)

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
# numSims = 100 #100
# numSeq = 50 #100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
# sigmasq = 0.25 #0.1

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
model0 = list(
  designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(
  designMat = desX1, beta.mean = mu1, beta.var = V1)

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenario 1: True function is quadratic
################################################################################
# betaT = c(-0.2, -0.4, 0.4)
betaT = c(0, 0, 1) # want quadratic coefficient to be much larger than linear
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
# curve(fT, from = xmin, to = xmax)

################################################################################
# run simulations
################################################################################

if(sequential_alpha | hybrid_alpha){
  seqmed_sims = foreach(i = 1:numSims) %dopar% {
    print(paste0("starting simulation ", i, " out of ", numSims))
    SeqMED(
      y.in = NULL, x.in = NULL, true.function = fT,
      model0 = model0, model1 = model1, 
      error.var = sigmasq, xmin = xmin, xmax = xmax,
      candidates = candidates, numSeq = numSeq, seqN = seqN, 
      alpha_seq = alpha_seq)
  }
} else{
  seqmed_sims = foreach(i = 1:numSims) %dopar% {
    print(paste0("starting simulation ", i, " out of ", numSims))
    SeqMED(
      y.in = NULL, x.in = NULL, true.function = fT,
      model0 = model0, model1 = model1, 
      error.var = sigmasq, xmin = xmin, xmax = xmax,
      candidates = candidates, numSeq = numSeq, seqN = seqN, 
      alpha_seq = alpha)
  }
}

seqmed_file0 = paste0(
  output_dir,
  "/sm", "_scen", scenario,
  "_N", Nttl,
  "_sigmasq", sigmasq,
  "_alpha", alpha
)
if(sequential_alpha){
  seqmed_file0 = paste0(seqmed_file0, "seq")
} else if(hybrid_alpha){
  seqmed_file0 = paste0(seqmed_file0, "hybrid")
}
seqmed_file = paste0(
  seqmed_file0, "_numSims", numSims, 
  ".rds")

saveRDS(seqmed_sims, file = seqmed_file)

