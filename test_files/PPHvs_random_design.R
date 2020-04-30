# testing higher dimensions

# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/add_MMEDvs.R", sep = ""))
source(paste(functions_home, "/SMMEDvs.R", sep = ""))

library(Matrix)
library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)

# --- evaluations --- #

source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))

mu_full = c(0.5, 0.5, 0.5) #

# settings
xmin = -1
xmax = 1
sigmasq01 = 0.5
sigmasq = 0.3
numCandidates = 1000 #
xmin = -1
xmax = 1
p = 3
k = 4 * p
initN = 5
pfull = length(mu_full)

alpha_seq = 1

# hypotheses
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]

# H0 true
betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]

# initial data
set.seed(12)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = 123))

numSeq_seq = c(50, 100, 250)
N_seq2 = 1
sigmasq2 = 0.3

models = list("H0" = list(mu0, V0, NULL, indices0), "H1" = list(mu1, V1, NULL, indices1))
par(mfrow=c(1,length(models)))
for(i in 1:length(numSeq_seq)){
  numSeq2 = numSeq_seq[i]
  Nnew2 = numSeq2 * N_seq2
  Ntot2 = initN + Nnew2
  # first get inputs
  x_random2 = matrix(runif(n = pfull * Nnew2, min = xmin, max = xmax), 
                     nrow = Nnew2, ncol = pfull)
  # then draw outputs from a multivariate normal
  y_random2 = simulateYvs(x_random2[ , indicesT], Nnew2, betaT, sigmasq2, 1, seed = 12)
  x_random2 = rbind(initD, x_random2)
  y_random2 = c(inity, y_random2)
  
  # calc EPPH
  EPPHseq_rand2 = calcEPPHseqdata(y_random2, x_random2, models, sigmasq2, initN, numSeq = numSeq2, N_seq = N_seq2)
  EPPHs_rand2 = calcEPPHdata(y_random2, x_random2, Ntot2, models, sigmasq2)
  
  for(k in 1:length(models)){
    if(k == 1){
      plot(x = 1:numSeq2, y = EPPHseq_rand2[k, ], type = "l", 
           xlab = "steps 1:numSeq", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
           ylim = c(0, 1), main = paste("n=", numSeq2, sep = ""))
      points(x = numSeq2, y = EPPHs_rand2[k])
    }
    if(k == 2){
      plot(x = 1:numSeq2, y = EPPHseq_rand2[k, ], type = "l", 
           xlab = "steps 1:numSeq", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
           ylim = c(0, 1))
      points(x = numSeq2, y = EPPHs_rand2[k])
    }
  }
}
