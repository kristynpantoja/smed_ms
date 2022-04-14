################################################################################
# last updated: 2/28/22
# purpose: to create a list of seqmed simulations
# scenario 1:
#   linear vs. quadratic,
#   where the true function is quadratic
# scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic
rm(list = ls())

scenario = 1 # 1, 2
given_Dinit = FALSE

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

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))

# for D-optimal design
library(AlgDesign)

# set up parallelization
rng.seed = 123 # 123, 345

library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggpubr)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100 # 1 simulation with 100 design points
numSeq = 100 # 100 design points
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
if(scenario == 1){
  if(Nttl == 100){
    sigmasq = 0.28
  } else if(Nttl == 12){
    sigmasq = 0.04
  }
} else if(scenario == 2){
  if(Nttl == 100){
    sigmasq = 0.21
  } else if(Nttl == 12){
    sigmasq = 0.038
  }
}
alpha = 1

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
# curve(fT, from = xmin, to = xmax)
W_fun = function(x, fT, error.var, beta.var){
  rho = error.var / beta.var
  fx = fT(x)
  squareddiff_means = (fx * 
                         ((1 + x^2) / (1 + x^2 + rho) - 
                            (1 + x^2 + x^4) / (1 + x^2 + x^4 + rho)))^2
  squareddiff_sds = error.var * (sqrt((1 + x^2) / (1 + x^2 + rho) + 1) -
                                sqrt((1 + x^2 + x^4) / (1 + x^2 + x^4 + rho)))^2
  return(squareddiff_means + squareddiff_sds)
}
W_fun_curve = function(x){
  W_fun(x = x, fT = fT, error.var = sigmasq, beta.var = sigmasq01)
}
curve(W_fun_curve, from = xmin, to = xmax)

