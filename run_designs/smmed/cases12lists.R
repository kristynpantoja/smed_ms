# libraries
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# source files for evaluations

# --- sources to generate MEDs --- #
home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED.R", sep = ""))

# --- sources to designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
source(paste(functions_home, "/plot_EPH1.R", sep = ""))
source(paste(functions_home, "/plot_MSE.R", sep = ""))
source(paste(functions_home, "/plot_posterior_variance.R", sep = ""))
source(paste(functions_home, "/postpredyhat_mse_closedform.R", sep = ""))

# CASE 1 #

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

# MED design #
typeT = 3
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
xmin = -1
xmax = 1

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
curve(f0, col = 2, lwd = 5, xlim = c(xmin, xmax), xlab = "x", ylab = "f")
curve(f1, add = T, col = 5, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)
legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))

# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^3
k = 4
xmin = -1
xmax = 1
p = 2
N = 100

numSeq = 10 
N_seq = 10
max_alpha = 2 * p
alpha_seq = c(0, ((2:numSeq) / numSeq) * (max_alpha))

numSeqMMED = 25
mmed_data_list = list()
for(i in 1:numSeqMMED){
  mmed_data_list[[i]] = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
                                        numCandidates, k, xmin, xmax, p, numSeq = numSeq, N_seq = N_seq, alpha_seq = alpha_seq, seed = 123 + i)
}
saveRDS(mmed_data_list, file = "case1smmeds.rds")


# CASE 2 #

# MED design #
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
betaT = c(0, -0.75, 0, 1)
typeT = 4
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
xmin = -1
xmax = 1
mu0 = c(0, 0)
V0 = diag(rep(sigmasq01,length(mu0)))
mu1 = c(0, 0, 0)
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
type01 = c(2, 3)
numCandidates = 10^3
k = 4
S = 10
xmin = -1
xmax = 1
p = 2
N = 100

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
curve(f0, col = 2, lwd = 5, xlim = c(xmin, xmax), xlab = "x", ylab = "f")
curve(f1, add = T, col = 5, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)
legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))

numSeq = 10 
N_seq = 10
max_alpha = 2 * p
alpha_seq = c(0, ((2:numSeq) / numSeq) * (max_alpha))

numSeqMMED = 25
mmed_data2_list = list()
for(i in 1:numSeqMMED){
  mmed_data2_list[[i]] = simulate_seqMED(betaT, typeT, mu0, mu1, V0, V1, sigmasq, f0, f1, type01, 
                                         numCandidates, k, xmin, xmax, p, numSeq = numSeq, 
                                         N_seq = N_seq, alpha_seq = alpha_seq, seed = 123 + i)
}
saveRDS(mmed_data2_list, file = "case2smmeds.rds")






