################################################################################
# last updated: 12/01/20
# purpose: to create a list of boxhill simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is quadratic

################################################################################
# Sources/Libraries
################################################################################
# tamu cluster
# home = "/scratch/user/kristynp/smed_ms/"
# output_home = paste(home,"run_designs/",sep="")
output_home = "run_designs/updated_simulations/lm"
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

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 25

# simulation settings
numSeq = 5
seqN = 5
N = numSeq * seqN
N.new = (numSeq - 1) * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)

# SeqMED settings
type01 = c(2, 3)
sigmasq = 0.1
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2

# boxhill settings
desX0 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x))
}
desX1 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x, x^2))
}

################################################################################
# Scenario 1: True function is quadratic
################################################################################
betaT = c(-0.2, -0.4, 0.4)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2

# boxhill settings
model0 = list(designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(designMat = desX1, beta.mean = mu1, beta.var = V1)
prior_probs = rep(1 / 2, 2)

# seqmed settings
typeT = 3

# generate boxhills
input_list = list()
bh_list = list()
for(i in 1:numSims){
  print(paste0("starting simulation ", i, " out of ", numSims))
  seqmed.res = SeqMED(
    D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
    mean_beta0 = mu0, mean_beta1 = mu1, var_beta0 = V0, var_beta1 = V1, 
    var_e = sigmasq, f0 = f0, f1 = f1, type = type01, 
    candidates = candidates, numSeq = 1, seqN = seqN, seed = 123 + i
    )
  x_input = seqmed.res$D
  y_input = seqmed.res$y
  input_list[[i]] = list(x_input = x_input, y_input = y_input)
  bh.res = BH_m2(y_input, x_input, prior_probs, model0, model1, N.new, 
                candidates, fT, sigmasq, seed = 1995 + i)
  bh_list[[i]] = bh.res
}
saveRDS(list(input_list = input_list, bh_list = bh_list), 
        paste(output_home, "/scenario1_boxhill_simulations", 
              "_N", N, 
              "_Nnew", N.new,
              ".rds", sep = ""))




