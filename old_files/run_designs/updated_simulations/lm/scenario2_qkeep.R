################################################################################
# last updated: 12/09/20
# purpose: to create a list of seqmed simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic

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
numSeq = 100
seqN = 1
N = numSeq * seqN
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
MMEDinputdata = FALSE
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
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenario 2: True function is cubic
################################################################################
betaT = c(0, -0.75, 0, 1)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3

# seqmed settings
typeT = 4

# generate seqmed
seqmed_newq = SeqMED(
  D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
  beta.mean0 = mu0, beta.mean1 = mu1, beta.var0 = V0, beta.var1 = V1, 
  error.var = sigmasq, f0 = f0, f1 = f1, type = type01, xmin = xmin, xmax = xmax, 
  candidates = candidates, numSeq = numSeq, seqN = seqN, seed = 123
)
seqmed_keepq = SeqMED(
  D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
  beta.mean0 = mu0, beta.mean1 = mu1, beta.var0 = V0, beta.var1 = V1, 
  error.var = sigmasq, f0 = f0, f1 = f1, type = type01, xmin = xmin, xmax = xmax, 
  candidates = candidates, numSeq = numSeq, seqN = seqN, seed = 123, 
  newq = FALSE
)
# all 100 points
plt1 = ggplot(data = data.frame(x = seqmed_newq$D), aes(x = x)) + 
  geom_histogram(binwidth = 0.08) + ggtitle("Original: re-calculate q")
plt2 = ggplot(data = data.frame(x = seqmed_keepq$D), aes(x = x)) + 
  geom_histogram(binwidth = 0.08) + ggtitle("Keep q")
ggarrange(plt1, plt2)

# first 3 points
plt1 = ggplot(data = data.frame(x = seqmed_newq$D[1:3]), aes(x = x)) + 
  geom_histogram(binwidth = 0.08) + ggtitle("Original: re-calculate q")
plt2 = ggplot(data = data.frame(x = seqmed_keepq$D[1:3]), aes(x = x)) + 
  geom_histogram(binwidth = 0.08) + ggtitle("Keep q")
ggarrange(plt1, plt2)

# first i points
i = 10
plt1 = ggplot(data = data.frame(x = seqmed_newq$D[1:i], 
                                y = 1:i * 0.1)) + 
  geom_histogram(binwidth = 0.08, aes(x = x), fill = "gray") +
  geom_text(aes(x = x, y = y, label = as.character(1:i)), color = 1) + 
  ggtitle("Original: re-calculate q") + 
  theme_classic() + ylim(0, 2)
plt2 = ggplot(data = data.frame(x = seqmed_keepq$D[1:i], 
                                y = 1:i * 0.1)) + 
  geom_histogram(binwidth = 0.08, aes(x = x), fill = "gray") +
  geom_text(aes(x = x, y = y, label = as.character(1:i)), color = 1) + 
  ggtitle("Keep q") + 
  theme_classic() + ylim(0, 2)
ggarrange(plt1, plt2)
