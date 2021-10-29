################################################################################
# last updated: 10/29/21
# purpose: to create a list of seqmed simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic

scenario = 2

beta_setting = 4 # 0, 1, 2, 3, 4

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lm/testing"
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
numSims = 100 #100
numSeq = 36 #36, 100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
sigmasq = 0.1 # 0.1

# shared settings
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
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
if(beta_setting == 0){
  betaT = c(0, -0.75, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 1){
  betaT = c(0, 0, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 2){
  betaT = c(0, 0, 1, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 3){
  betaT = c(0, 0, 1)
  fx_outer = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
  discontinuity = 0.01
  height = 1
  fx_innerneg = function(x) height - ((fx_outer(discontinuity) - height) / discontinuity) * x
  fx_innerpos = function(x) height + ((fx_outer(discontinuity) - height) / discontinuity) * x
  fT = function(x){
    fx = rep(NA, length(x))
    for(i in 1:length(x)){
      if(abs(x[i]) >= discontinuity){
        fx[i] = fx_outer(x[i])
      } else if(x[i] > -discontinuity & x[i] <= 0){
        fx[i] = fx_innerneg(x[i])
      } else{
        fx[i] = fx_innerpos(x[i])
      }
    }
    return(fx)
  }
} else if(beta_setting == 4){
  betaT = c(0, 0, 1)
  fx_outer = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
  discontinuity = 0.01
  height = 1
  fx_inner = function(x) height
  fT = function(x){
    fx = rep(NA, length(x))
    for(i in 1:length(x)){
      if(abs(x[i]) >= discontinuity){
        fx[i] = fx_outer(x[i])
      } else{
        fx[i] = fx_inner(x[i])
      }
    }
    return(fx)
  }
}
curve(fT, from = xmin, to = xmax)

################################################################################
# run simulations
################################################################################

# generate seqmeds
seqmed_sims = foreach(i = 1:numSims) %dopar% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  # SeqMED(
  #   D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
  #   beta.mean0 = mu0, beta.mean1 = mu1, beta.var0 = V0, beta.var1 = V1, 
  #   error.var = sigmasq, f0 = f0, f1 = f1, type = type01, xmin = xmin, xmax = xmax, 
  #   candidates = candidates, numSeq = numSeq, seqN = seqN)
  SeqMED(
    y.in = NULL, x.in = NULL, true.function = fT,
    model0 = model0, model1 = model1, 
    error.var = sigmasq, xmin = xmin, xmax = xmax,
    candidates = candidates, numSeq = numSeq, seqN = seqN)
}

saveRDS(seqmed_sims, file = paste0(
  output_dir,
  "/sm", "_scen", scenario, 
  "_sigmasq", sigmasq,
  "_Nttl", Nttl, 
  "_numSims", numSims, 
  "_beta", beta_setting,
  ".rds"
))







