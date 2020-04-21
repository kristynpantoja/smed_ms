# --- Working Directory --- #

# Computer
# home = "/home/kristyn/Documents/smed_ms"
# output_home = paste(home, "/", sep = "")

# Cluster
home = "/scratch/user/kristynp/smed_ms"
output_home = paste(home,"/run_designs_lmvs7/",sep="")

jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/generate_MMEDgp_oneatatime.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)

# --- sources to generate MEDs --- #
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED_multidim.R", sep = ""))

# --- simulations  --- #
numSims = 1

# --- settings  --- #

mu_full = c(0.5, 0.5, 0.5, 0.5, 0.5) #

# settings
xmin = -1
xmax = 1
sigmasq01 = 0.25
sigmasq = 0.3
numCandidates = 10000 #
xmin = -1
xmax = 1
p = 3
k = 4 * p
initN = 50
pfull = length(mu_full)

# sequential settings
numSeq = 10
N_seq = 10 #
alpha_seq = 1

# hypotheses
indices0 = c(1, 2, 3) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]

# sim for H0 case
betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]
seed = jid + 1
set.seed(seed)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateY_multidim(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
smmed_H0 = simulate_seqMED_multidim(mu_full, betaT, indicesT, indices0, indices1, mu0, mu1, 
                                             sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
                                             initD, inity, numSeq, N_seq, seed = seed)
saveRDS(smmed_H0, paste(output_home, "lmvs3v5H0_sim", jid, ".rds", sep = ""))

# sim for H1 case
betaT = mu_full[indices1]
indicesT = indices1
fT = function(x) betaT %*% x[indicesT]
seed = jid + 1
set.seed(seed)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateY_multidim(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
smmed_H1 = simulate_seqMED_multidim(mu_full, betaT, indicesT, indices0, indices1, mu0, mu1, 
                                             sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
                                             initD, inity, numSeq, N_seq, seed = seed)
saveRDS(smmed_H1, paste(output_home, "lmvs3v5H1_sim", jid, ".rds", sep = ""))
