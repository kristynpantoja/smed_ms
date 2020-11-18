################################################################################
### --- Linear Models ------------------------------------------------------ ###
################################################################################

# --- Working Directory --- #
home = "/home/kristyn/Documents/research/seqmed/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
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

# for plots
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
image_path = "/home/kristyn/Pictures"

### --- shared settings for both scenarios : linear vs. quadratic --- ###

# --- simulations --- #

# settings: Squared Exponential vs. Matern
N = 10
N.new = 40
numSeq = 5; N_seq = 10
xmin = -1
xmax = 1
numCandidates = 10^3
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)

type01 = c(2, 3)
sigmasq = 0.1
betaT = c(-0.2, -0.4, 0.4)
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))

f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2

################################################################################
# input data
################################################################################

# input points : randomly-selected points
set.seed(1997)
x_input = runif(N, xmin, xmax)

################################################################################
# Scenario 1: True function is quadratic
################################################################################

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
typeT = 3
y_input = fT(x_input) + rnorm(N, 0, sqrt(sigmasq))

# seqmed
seqmed.res = SeqMED(
  D1 = x_input, y1 = y_input, true_beta = betaT, true_type = typeT, 
  mean_beta0 = mu0, mean_beta1 = mu1, var_beta0 = V0, var_beta1 = V1, 
  var_e = sigmasq, f0 = f0, f1 = f1, type = type01, 
  candidates = candidates, numSeq = numSeq, N_seq = N_seq
)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(seqmed.res$D, seqmed.res$y, col = 2, cex = 0.5)

# boxhill
prior_probs = rep(1 / 2, 2)
model0 = list(
  X = constructDesignX(x_input, N, type = type01[1]), 
  X.n = constructDesignX(candidates, numCandidates, type = type01[1]),
  beta.mean = mu0, 
  beta.var = V0
)
model1 = list(
  X = constructDesignX(x_input, N, type = type01[2]), 
  X.n = constructDesignX(candidates, numCandidates, type = type01[2]),
  beta.mean = mu1, 
  beta.var = V1
)
BH_m2(y = y_input, prior.probs = prior_probs, error.var = sigmasq, 
      model.i = model0, model.j = model1)
  




