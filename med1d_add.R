library(mvtnorm)

# --- sources to generate MEDs --- #
home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))

source(paste(functions_home, "/simulate_y.R", sep = ""))

# --- example 2.1 --- #

# MED design #
sigmasq01 = 0.01
mu0 = c(0, 0)
V0 = diag(rep(sigmasq01,length(mu0)))
mu1 = c(-0.2, 0, 0.5)
V1 = diag(rep(sigmasq01,length(mu1)))
sigmasq = 0.01 # 1
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type = c(2, 3)
numCandidates = 10^3
k = 4
S = 10
xmin = -1
xmax = 1
p = 2
Ninit = 50
# alpha = 2p
# design_ex2.1.1_k4_alpha2p = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type, Ninit,
#                               numCandidates, k, xmin, xmax, p = p, alpha = NULL)
# saveRDS(design_ex2.1.1_k4_alpha2p, file = "med1d_ex2pt1pt1_k4_alpha2p.rds")
design_MED_alpha2p = readRDS("med1d_ex2pt1pt1_k4_alpha2p.rds")
hist(design_MED_alpha2p, breaks = 20)
# original MED
# design_ex2.1.1_k4_alpha1 = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type, Ninit,
#                               numCandidates, k, xmin, xmax, p = p, alpha = 1)
# saveRDS(design_ex2.1.1_k4_alpha1, file = "med1d_ex2pt1pt1_k4_alpha1.rds")
design_MED_alpha1 = readRDS("med1d_ex2pt1pt1_k4_alpha1.rds")
hist(design_MED_alpha1, breaks = 20)

# add points (without data)
design_ex2.1.1_k4_alpha2p_update = add_MED_ms_oneatatime(design_MED_alpha2p, mu0, mu1, V0, V1, sigmasq,
                                                         f0, f1, type, N2 = 50, numCandidates, k, 
                                                         xmin, xmax, p, alpha = NULL)
hist(design_ex2.1.1_k4_alpha2p_update$updatedD, breaks = 20)
design_MED_alpha2p_nodata = design_ex2.1.1_k4_alpha2p_update$updatedD
# see if it's the same as if the 100 points were generated all at once
# Nttl = 100
# design_ex2.1.1_k4_alpha2p = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type, Nttl,
#                                               numCandidates, k, xmin, xmax, p = p, alpha = NULL)
# hist(design_ex2.1.1_k4_alpha2p, breaks = 20)
# all.equal(design_ex2.1.1_k4_alpha2p_update$updatedD, design_ex2.1.1_k4_alpha2p) # TRUE - yay!

# add points (with data)
betaT = mu1
typeT = 3
y = as.vector(simulateY(design_MED_alpha2p, N, mu1, sigmasq, 1, typeT, seed = 123))
plot(x = design_MED_alpha2p, y = y)

design_ex2.1.1_k4_alpha2p_update = add_MED_ms_oneatatime_data(design_MED_alpha2p, y, mu0, mu1, V0, V1, sigmasq,
                                                         f0, f1, type, N2 = 50, numCandidates, k, 
                                                         xmin, xmax, p, alpha = NULL)
hist(design_ex2.1.1_k4_alpha2p_update$updatedD, breaks = 20)
hist(design_MED_alpha2p_nodata, breaks = 20)

length(design_ex2.1.1_k4_alpha2p_update$q_initD)




