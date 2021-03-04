# --- Working Directory --- #

# Computer
# home = "/home/kristyn/Documents/smed_ms"
# output_home = paste(home, "/run_designs/smmed/run_designs_lmvs/", sep = "")
# jid = 5

# Cluster
home = "/scratch/user/kristynp/smed_ms"
output_home = paste(home,"/run_designs_lmvs/",sep="")
jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

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

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)

# --- simulations  --- #
numSims = 1

# --- settings  --- #

mu_full = c(0.5, 0.5, 0.5) #

# settings
xmin = -1
xmax = 1
sigmasq01 = 0.5
sigmasq = 0.3
numCandidates = 5000 #
xmin = -1
xmax = 1
p = 3
k = 4 * p
initN = 5
pfull = length(mu_full)

# sequential settings
numSeq = 9
N_seq = 3
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

# sim for H0 case
betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]
seed = jid + 1
set.seed(seed)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
smmed_H0 = generate_SMMEDvs(mu_full, betaT, indicesT, indices0, indices1, mu0, mu1,
                            sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
                            initD, inity, numSeq, N_seq, seed = 12)
saveRDS(smmed_H0, paste(output_home, "lmvs2v3H0_sim", jid, ".rds", sep = ""))

# sim for H1 case
betaT = mu_full[indices1]
indicesT = indices1
fT = function(x) betaT %*% x[indicesT]
seed = jid + 1
set.seed(seed)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = seed))
smmed_H1 = generate_SMMEDvs(mu_full, betaT, indicesT, indices0, indices1, mu0, mu1,
                            sigmasq, sigmasq01, V0, V1, xmin, xmax, numCandidates, k, p,
                            initD, inity, numSeq, N_seq, seed = 12)
saveRDS(smmed_H1, paste(output_home, "lmvs2v3H1_sim", jid, ".rds", sep = ""))
