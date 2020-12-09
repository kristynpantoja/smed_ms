################################################################################
### --- Linear Models ------------------------------------------------------ ###
################################################################################

# --- Working Directory --- #
# home = "/home/kristyn/Documents/research/seqmed/smed_ms"

# --- Sources/Libraries --- #
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
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)

type01 = c(2, 3)
sigmasq = 0.1
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

f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2

################################################################################
# input data
################################################################################

# input points : randomly-selected points
set.seed(1997)
x_input = runif(N, xmin, xmax)
# x_input = c(rep(-1, N / 2), rep(1, N / 2))
# x_input = rep(0, N)

################################################################################
# Scenario 1: True function is quadratic
################################################################################
betaT = c(-0.2, -0.4, 0.4)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
typeT = 3
y_input = fT(x_input) + rnorm(N, 0, sqrt(sigmasq))

################################################################################
# seqmed
seqmed.res = SeqMED(
  D1 = x_input, y1 = y_input, true_beta = betaT, true_type = typeT, 
  mean_beta0 = mu0, mean_beta1 = mu1, var_beta0 = V0, var_beta1 = V1, 
  var_e = sigmasq, f0 = f0, f1 = f1, type = type01, 
  candidates = candidates, numSeq = numSeq, seqN = N_seq, 
  seed = 1
)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(seqmed.res$D, seqmed.res$y, col = 2, cex = 0.5)

# seqmed without input data
seqmed.res = SeqMED(
  D1 = NULL, y1 = NULL, true_beta = betaT, true_type = typeT, 
  mean_beta0 = mu0, mean_beta1 = mu1, var_beta0 = V0, var_beta1 = V1, 
  var_e = sigmasq, f0 = f0, f1 = f1, type = type01, 
  candidates = candidates, numSeq = numSeq, seqN = N_seq, 
  seed = 1
)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(seqmed.res$D, seqmed.res$y, col = 2, cex = 0.5)

################################################################################
# box and hill method
model0 = list(designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(designMat = desX1, beta.mean = mu1, beta.var = V1)
prior_probs = rep(1 / 2, 2)

BHres = BH_m2(y_input, x_input, prior_probs, model0, model1, N.new, candidates, 
              fT, sigmasq, 
              seed = 1)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(BHres$x.new, BHres$y.new, col = 2, cex = 0.5)
BHres$x.new

length(which(BHres$x.new < -0.5))
length(which(BHres$x.new > 0.5))
length(which(BHres$x.new >= -0.5 & BHres$x.new <= 0.5))

# box hill without preliminary data
BHres = BH_m2(NULL, NULL, prior_probs, model0, model1, N + N.new, candidates, 
              fT, sigmasq, 
              seed = 1)
plot(BHres$x, BHres$y, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(BHres$x.new, BHres$y.new, col = 2, cex = 0.5)

################################################################################
# Scenario 2: True function is cubic
################################################################################
betaT = c(0, -0.75, 0, 1)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
typeT = 4
y_input = fT(x_input) + rnorm(N, 0, sqrt(sigmasq))

################################################################################
# seqmed
seqmed.res = SeqMED(
  D1 = x_input, y1 = y_input, true_beta = betaT, true_type = typeT, 
  mean_beta0 = mu0, mean_beta1 = mu1, var_beta0 = V0, var_beta1 = V1, 
  var_e = sigmasq, f0 = f0, f1 = f1, type = type01, 
  candidates = candidates, numSeq = numSeq, seqN = N_seq, 
  seed = 1
)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(seqmed.res$D, seqmed.res$y, col = 2, cex = 0.5)

################################################################################
# box and hill method
model0 = list(designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(designMat = desX1, beta.mean = mu1, beta.var = V1)

# calculate prior probabilities using preliminary data (input data)
prior_probs = rep(1 / 2, 2)

BHres = BH_m2(y_input, x_input, prior_probs, model0, model1, N.new, candidates, 
              fT, sigmasq, 
              seed = 1)
plot(x_input, y_input, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(BHres$x.new, BHres$y.new, col = 2, cex = 0.5)
BHres$x.new

length(which(BHres$x.new < -0.5))
length(which(BHres$x.new > 0.5))
length(which(BHres$x.new >= -0.5 & BHres$x.new <= 0.5))

# box hill without preliminary data
BHres = BH_m2(NULL, NULL, prior_probs, model0, model1, N + N.new, candidates, 
              fT, sigmasq, 
              seed = 1)
plot(BHres$x, BHres$y, ylim = c(-2, 2), xlim = c(-2, 2))
curve(fT, add = TRUE)
points(BHres$x.new, BHres$y.new, col = 2, cex = 0.5)

