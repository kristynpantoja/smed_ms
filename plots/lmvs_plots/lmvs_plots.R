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
source(paste(functions_home, "/posterior_mean_mse.R", sep = ""))


  
# --- settings for 2 vs 3 factors --- #
### These are shared for both H1 true and H2 true cases ###

mu_full = c(0.5, 0.5, 0.5) #

# settings
xmin = -1
xmax = 1
sigmasq01 = 0.5
sigmasq = 0.5
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

###################
# --- H0 true --- #
###################

betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]
set.seed(12)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = 123))

# import H0 sims
numSimsH0 = 50
simsH0 = list()
for(i in 1:numSimsH0){
  simsH0[[i]] = readRDS(paste(home, "/run_designs/smmed/run_designs_lmvs3/lmvs2v3H0_sim", i, ".rds", sep = ""))
}
smmedvsH0 = simsH0[[1]]

Nnew = numSeq * N_seq
Ntot = Nnew + initN

# --- other designs --- #

# random design
set.seed(12345)
x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                  nrow = Nnew, ncol = pfull)
y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_random = rbind(initD, x_random)
y_random = c(inity, y_random)

# doptimal design
dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
for(i in 0:(Nnew - 1)){
  x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
}
y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_doptimal = rbind(initD, x_doptimal)
y_doptimal = c(inity, y_doptimal)

# 3 level factorial design
factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
###
x_factorial = factorial_pts
###
y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_factorial = rbind(initD, x_factorial)
y_factorial = c(inity, y_factorial)



# --- EPPH non-sequential --- #

numSims_pph = 50
models = list("H0" = list(mu0, V0, NULL, indices0),
              "H1" = list(mu1, V1, NULL, indices1))
#EPPHs for other designs
EPPHs_rand = calcEPPH(x_random, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_dopt = calcEPPH(x_doptimal, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_fact = calcEPPH(x_factorial, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
#EPPHs for smmed
EPPHs_smmed = calcEPPHseqdata(smmedvsH0$y, smmedvsH0$D, models, sigmasq, initN, numSeq, N_seq)

#EPPHs for smmed sims
idxlast = numSeq + 1
EPPH0_smmedsims = matrix(NA, idxlast, numSimsH0)
EPPH1_smmedsims = matrix(NA, idxlast, numSimsH0)
for(i in 1:numSimsH0){
  simH0.temp = simsH0[[i]]
  EPPH_temp = calcEPPHseqdata(simH0.temp$y, simH0.temp$D, models, sigmasq, initN, numSeq, N_seq)
  EPPH0_smmedsims[ , i] = EPPH_temp[1 , ]
  EPPH1_smmedsims[ , i] = EPPH_temp[2 , ]
}

EPPH0_smmedsims_mean = apply(EPPH0_smmedsims, 1, mean)
EPPH1_smmedsims_mean = apply(EPPH1_smmedsims, 1, mean)

EPPH0_smmedsims_med = apply(EPPH0_smmedsims, 1, median)
EPPH1_smmedsims_med = apply(EPPH1_smmedsims, 1, median)

design_names = c("smmed", "random", "dopt", "factorial")
design_col = c(1, 3, 4, 6)

par(mfrow=c(1,length(models)))
#mean
for(k in 1:length(models)){
  if(k == 1){
    plot(x = 1:idxlast, y = EPPH0_smmedsims_mean, type = "l", 
         xlab = "", ylab = paste("E[P(H", k - 1, "|Y,X)|X]", sep = ""), 
         ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
    points(x = idxlast, y = EPPH0_smmedsims_mean[idxlast])
  }
  if(k == 2){
    plot(x = 1:idxlast, y = EPPH1_smmedsims_mean, type = "l", 
         xlab = "", ylab = paste("E[P(H", k - 1, "|Y,X)|X]", sep = ""), 
         ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
    points(x = idxlast, y = EPPH1_smmedsims_mean[idxlast])
  }
  abline(h = EPPHs_rand[k], col = 3, lty = 3)
  abline(h = EPPHs_dopt[k], col = 4, lty = 3)
  abline(h = EPPHs_fact[k], col = 6, lty = 3)
  points(x = idxlast, y = EPPHs_rand[k], col = 3, pch = 16)
  points(x = idxlast, y = EPPHs_dopt[k], col = 4, pch = 17)
  points(x = idxlast, y = EPPHs_fact[k], col = 6, pch = 18)
  if(k == 1) legend("bottomright", design_names, lty = c(1, rep(3, 3)), pch = c(1, 16:18), col = design_col, bg = "white")
}
# # median
# for(k in 1:length(models)){
#   if(k == 1){
#     plot(x = 1:numSeq, y = EPPH0_smmedsims_med, type = "l", 
#          xlab = "# steps", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
#          ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
#     points(x = numSeq, y = EPPH0_smmedsims_med[numSeq])
#   }
#   if(k == 2){
#     plot(x = 1:numSeq, y = EPPH1_smmedsims_med, type = "l", 
#          xlab = "# steps", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
#          ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
#     points(x = numSeq, y = EPPH1_smmedsims_med[numSeq])
#   }
#   abline(h = EPPHs_rand[k], col = 3, lty = 3)
#   abline(h = EPPHs_dopt[k], col = 4, lty = 3)
#   abline(h = EPPHs_fact[k], col = 6, lty = 3)
#   points(x = numSeq, y = EPPHs_rand[k], col = 3, pch = 16)
#   points(x = numSeq, y = EPPHs_dopt[k], col = 4, pch = 17)
#   points(x = numSeq, y = EPPHs_fact[k], col = 6, pch = 18)
#   if(k == 1) legend("bottomright", design_names, lty = c(1, rep(3, 3)), pch = c(1, 16:18), col = design_col, bg = "white")
# }



# --- EPPH sequential --- #

numSims_seqPPH = numSimsH0

rand_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
dopt_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
fact_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
smmed_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)

rand_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
dopt_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
fact_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
smmed_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)

for(j in 1:numSims_seqPPH){
  set.seed(j + 50)
  smmedsim = simsH0[[j]]
  initD = smmedsim$initD
  inity = smmedsim$y[1:initN]
  
  # random design
  x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                    nrow = Nnew, ncol = pfull)
  y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_random = rbind(initD, x_random)
  y_random = c(inity, y_random)
  
  # doptimal design
  dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
  x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
  for(i in 0:(Nnew - 1)){
    x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
  }
  y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_doptimal = rbind(initD, x_doptimal)
  y_doptimal = c(inity, y_doptimal)
  
  # 3 level factorial design
  factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
  factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
  ###
  x_factorial = factorial_pts
  ###
  y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_factorial = rbind(initD, x_factorial)
  y_factorial = c(inity, y_factorial)
  
  rand_seqPPH.temp = calcEPPHseqdata(y_random, x_random, models, sigmasq, initN, numSeq, N_seq)
  rand_seqPPH0[ , j] = rand_seqPPH.temp[1, ]
  rand_seqPPH1[ , j] = rand_seqPPH.temp[2, ]
  dopt_seqPPH.temp = calcEPPHseqdata(y_doptimal, x_doptimal, models, sigmasq, initN, numSeq, N_seq)
  dopt_seqPPH0[ , j] = dopt_seqPPH.temp[1, ]
  dopt_seqPPH1[ , j] = dopt_seqPPH.temp[2, ]
  fact_seqPPH.temp = calcEPPHseqdata(y_factorial, x_factorial, models, sigmasq, initN, numSeq, N_seq)
  fact_seqPPH0[ , j] = fact_seqPPH.temp[1, ]
  fact_seqPPH1[ , j] = fact_seqPPH.temp[2, ]
  smmed_seqPPH.temp = calcEPPHseqdata(smmedsim$y, smmedsim$D, models, sigmasq, initN, numSeq, N_seq)
  smmed_seqPPH0[ , j] = smmed_seqPPH.temp[1, ]
  smmed_seqPPH1[ , j] = smmed_seqPPH.temp[2, ]
}

par(mfrow=c(1,2))

#mean
rand_seqPPH0_mean = apply(rand_seqPPH0, 1, mean)
dopt_seqPPH0_mean = apply(dopt_seqPPH0, 1, mean)
fact_seqPPH0_mean = apply(fact_seqPPH0, 1, mean)
smmed_seqPPH0_mean = apply(smmed_seqPPH0, 1, mean)
rand_seqPPH1_mean = apply(rand_seqPPH1, 1, mean)
dopt_seqPPH1_mean = apply(dopt_seqPPH1, 1, mean)
fact_seqPPH1_mean = apply(fact_seqPPH1, 1, mean)
smmed_seqPPH1_mean = apply(smmed_seqPPH1, 1, mean)
# PPH0
idxlast = numSeq + 1
plot(x = 1:idxlast, y = smmed_seqPPH0_mean, type = "l",
     xlab = "", ylab = "E[P(H0|Y,X)|X]", ylim = c(0, 1))
lines(x = 1:idxlast, rand_seqPPH0_mean, col = 3, lty = 2)
lines(x = 1:idxlast, dopt_seqPPH0_mean, col = 4, lty = 2)
lines(x = 1:idxlast, fact_seqPPH0_mean, col = 6, lty = 2)
points(x = idxlast, y = rand_seqPPH0_mean[idxlast], col = 3, pch = 16)
points(x = idxlast, y = dopt_seqPPH0_mean[idxlast], col = 4, pch = 17)
points(x = idxlast, y = fact_seqPPH0_mean[idxlast], col = 6, pch = 18)
points(x = idxlast, y = smmed_seqPPH0_mean[idxlast])
legend("bottomright", design_names, lty = c(1, rep(2, length(design_names) - 1)), 
       pch = c(1, 16, 17, 18), col = design_col, bg = "white")
# PPH1
plot(x = 1:idxlast, y = smmed_seqPPH1_mean, type = "l",
     xlab = "", ylab = "E[P(H1|Y,X)|X]", ylim = c(0, 1))
lines(x = 1:idxlast, rand_seqPPH1_mean, col = 3, lty = 2)
lines(x = 1:idxlast, dopt_seqPPH1_mean, col = 4, lty = 2)
lines(x = 1:idxlast, fact_seqPPH1_mean, col = 6, lty = 2)
points(x = idxlast, y = rand_seqPPH1_mean[idxlast], col = 3, pch = 16)
points(x = idxlast, y = dopt_seqPPH1_mean[idxlast], col = 4, pch = 17)
points(x = idxlast, y = fact_seqPPH1_mean[idxlast], col = 6, pch = 18)
points(x = idxlast, y = smmed_seqPPH1_mean[idxlast])

# 
# # median
# rand_seqPPH0_median = apply(rand_seqPPH0, 1, median)
# dopt_seqPPH0_median = apply(dopt_seqPPH0, 1, median)
# fact_seqPPH0_median = apply(fact_seqPPH0, 1, median)
# smmed_seqPPH0_median = apply(smmed_seqPPH0, 1, median)
# rand_seqPPH1_median = apply(rand_seqPPH1, 1, median)
# dopt_seqPPH1_median = apply(dopt_seqPPH1, 1, median)
# fact_seqPPH1_median = apply(fact_seqPPH1, 1, median)
# smmed_seqPPH1_median = apply(smmed_seqPPH1, 1, median)
# # PPH0
# plot(x = 1:numSeq, y = smmed_seqPPH0_median, type = "l",
#      xlab = "steps 1:numSeq", ylab = "P(H0|Y)", ylim = c(0, 1), 
#      main = "Median P(H0|Y)")
# lines(x = 1:numSeq, rand_seqPPH0_median, col = 3, lty = 2)
# lines(x = 1:numSeq, dopt_seqPPH0_median, col = 4, lty = 2)
# lines(x = 1:numSeq, fact_seqPPH0_median, col = 6, lty = 2)
# points(x = numSeq, y = rand_seqPPH0_median[numSeq], col = 3, pch = 16)
# points(x = numSeq, y = dopt_seqPPH0_median[numSeq], col = 4, pch = 17)
# points(x = numSeq, y = fact_seqPPH0_median[numSeq], col = 6, pch = 18)
# points(x = numSeq, y = smmed_seqPPH0_median[numSeq])
# # PPH1
# plot(x = 1:numSeq, y = smmed_seqPPH1_median, type = "l",
#      xlab = "steps 1:numSeq", ylab = "P(H1|Y)", ylim = c(0, 1), 
#      main = "Median P(H1|Y)")
# lines(x = 1:numSeq, rand_seqPPH1_median, col = 3, lty = 2)
# lines(x = 1:numSeq, dopt_seqPPH1_median, col = 4, lty = 2)
# lines(x = 1:numSeq, fact_seqPPH1_median, col = 6, lty = 2)
# points(x = numSeq, y = rand_seqPPH1_median[numSeq], col = 3, pch = 16)
# points(x = numSeq, y = dopt_seqPPH1_median[numSeq], col = 4, pch = 17)
# points(x = numSeq, y = fact_seqPPH1_median[numSeq], col = 6, pch = 18)
# points(x = numSeq, y = smmed_seqPPH1_median[numSeq])



# --- MSE(Bn) --- #

par(mfrow = c(1, length(betaT)))

smmedvs = smmedvsH0
hyp_mu = mu0
hyp_V = V0
hyp_ind = indices0
mseBn_smmed = getClosedMSE(smmedvs$D, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_rand = getClosedMSE(x_random, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_dopt = getClosedMSE(x_doptimal, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_fact = getClosedMSE(x_factorial, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean

for(i in 1:length(betaT)){
  barplot(c(mseBn_smmed[i], mseBn_rand[i], mseBn_dopt[i], mseBn_fact[i]), names.arg = design_names, las = 2)
}


# --- the design --- #

par(mfrow = c(1, 3))

numBreaks = 8
smmed = smmedvsH0
Ntot = dim(smmed$D)[1]

maxcounts = rep(NA, length(mu_full))
for(i in 1:length(mu_full)){
  marginal = i
  h = hist(smmed$D[ (initN + 1):Ntot, marginal], breaks = numBreaks, plot = FALSE)
  maxcounts[i] = max(h$counts)
}

for(i in 1:length(mu_full)){
  marginal = i
  h1 = hist(smmed$D[ (initN + 1):Ntot, marginal], breaks = numBreaks, 
            ylim = range(0, maxcounts), main = "", 
            xlab = paste("marginal ", i, sep = ""))
}


###################
# --- H1 true --- #
###################

betaT = mu_full[indices1]
indicesT = indices1
fT = function(x) betaT %*% x[indicesT]
set.seed(12)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = 123))

# import H1 sims
numSimsH1 = 50
simsH1 = list()
for(i in 1:numSimsH1){
  simsH1[[i]] = readRDS(paste(home, "/run_designs/smmed/run_designs_lmvs3/lmvs2v3H1_sim", i, ".rds", sep = ""))
}
smmedvsH1 = simsH1[[1]]

# --- other designs --- #

Nnew = numSeq * N_seq
Ntot = Nnew + initN

# random design
set.seed(123)
x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                  nrow = Nnew, ncol = pfull)
y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_random = rbind(initD, x_random)
y_random = c(inity, y_random)

# postmean1 of random design
postmean0_random = matrix(NA, nrow = length(mu0), numSeq)
postmean1_random = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_random[ , i] = postmean(y_random[ind_seq], x_random[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_random[ , i] = postmean(y_random[ind_seq], x_random[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
  
}

# doptimal design
dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
for(i in 0:(Nnew - 1)){
  x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
}
y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_doptimal = rbind(initD, x_doptimal)
y_doptimal = c(inity, y_doptimal)

# postmean1 of doptimal design
postmean0_doptimal = matrix(NA, nrow = length(mu0), numSeq)
postmean1_doptimal = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_doptimal[ , i] = postmean(y_doptimal[ind_seq], x_doptimal[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_doptimal[ , i] = postmean(y_doptimal[ind_seq], x_doptimal[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
}

# 3 level factorial design
factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
###
x_factorial = factorial_pts
###
y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_factorial = rbind(initD, x_factorial)
y_factorial = c(inity, y_factorial)

# postmean1 of doptimal design
postmean0_factorial = matrix(NA, nrow = length(mu0), numSeq)
postmean1_factorial = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_factorial[ , i] = postmean(y_factorial[ind_seq], x_factorial[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_factorial[ , i] = postmean(y_factorial[ind_seq], x_factorial[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
}



# -- EPPH non-sequential --- #

numSims_pph = 50
models = list("H0" = list(mu0, V0, NULL, indices0),
              "H1" = list(mu1, V1, NULL, indices1))
#EPPHs for other designs
EPPHs_rand = calcEPPH(x_random, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_dopt = calcEPPH(x_doptimal, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_fact = calcEPPH(x_factorial, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
#EPPHs for smmed
EPPHs_smmed = calcEPPHseqdata(smmedvsH1$y, smmedvsH1$D, models, sigmasq, initN, numSeq, N_seq)

#EPPHs for smmed sims
idxlast = numSeq + 1
EPPH0_smmedsims = matrix(NA, idxlast, numSimsH1)
EPPH1_smmedsims = matrix(NA, idxlast, numSimsH1)
for(i in 1:numSimsH1){
  simH1.temp = simsH1[[i]]
  EPPH_temp = calcEPPHseqdata(simH1.temp$y, simH1.temp$D, models, sigmasq, initN, numSeq, N_seq)
  EPPH0_smmedsims[ , i] = EPPH_temp[1 , ]
  EPPH1_smmedsims[ , i] = EPPH_temp[2 , ]
}
EPPH0_smmedsims_mean = apply(EPPH0_smmedsims, 1, mean)
EPPH1_smmedsims_mean = apply(EPPH1_smmedsims, 1, mean)

EPPH0_smmedsims_med = apply(EPPH0_smmedsims, 1, median)
EPPH1_smmedsims_med = apply(EPPH1_smmedsims, 1, median)

design_names = c("smmed", "random", "dopt", "factorial")
design_col = c(1, 3, 4, 6)

par(mfrow=c(1,length(models)))
#mean
for(k in 1:length(models)){
  if(k == 1){
    plot(x = 1:idxlast, y = EPPH0_smmedsims_mean, type = "l", 
         xlab = "", ylab = paste("E[P(H", k - 1, "|Y,X)|X]", sep = ""), 
         ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
    points(x = idxlast, y = EPPH0_smmedsims_mean[idxlast])
  }
  if(k == 2){
    plot(x = 1:idxlast, y = EPPH1_smmedsims_mean, type = "l", 
         xlab = "", ylab = paste("E[P(H", k - 1, "|Y,X)|X]", sep = ""), 
         ylim = range(EPPHs_smmed[ k, ], EPPHs_rand[k], EPPHs_dopt[k], 0, 1))
    points(x = idxlast, y = EPPH1_smmedsims_mean[idxlast])
  }
  abline(h = EPPHs_rand[k], col = 3, lty = 3)
  abline(h = EPPHs_dopt[k], col = 4, lty = 3)
  abline(h = EPPHs_fact[k], col = 6, lty = 3)
  points(x = idxlast, y = EPPHs_rand[k], col = 3, pch = 16)
  points(x = idxlast, y = EPPHs_dopt[k], col = 4, pch = 17)
  points(x = idxlast, y = EPPHs_fact[k], col = 6, pch = 18)
  if(k == 1) legend("topright", design_names, lty = c(1, rep(3, 3)), pch = c(1, 16:18), col = design_col, bg = "white")
}



# --- EPPH sequential --- #

numSims_seqPPH = numSimsH1

idxlast = numSeq + 1
rand_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
dopt_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
fact_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
smmed_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)

rand_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
dopt_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
fact_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
smmed_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)

for(j in 1:numSims_seqPPH){
  set.seed(j + 50)
  smmedsim = simsH1[[j]]
  initD = smmedsim$initD
  inity = smmedsim$y[1:initN]
  
  # random design
  x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                    nrow = Nnew, ncol = pfull)
  y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_random = rbind(initD, x_random)
  y_random = c(inity, y_random)
  
  # doptimal design
  dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
  x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
  for(i in 0:(Nnew - 1)){
    x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
  }
  y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_doptimal = rbind(initD, x_doptimal)
  y_doptimal = c(inity, y_doptimal)
  
  # 3 level factorial design
  factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
  factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
  ###
  x_factorial = factorial_pts
  ###
  y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_factorial = rbind(initD, x_factorial)
  y_factorial = c(inity, y_factorial)
  
  rand_seqPPH.temp = calcEPPHseqdata(y_random, x_random, models, sigmasq, initN, numSeq, N_seq)
  rand_seqPPH0[ , j] = rand_seqPPH.temp[1, ]
  rand_seqPPH1[ , j] = rand_seqPPH.temp[2, ]
  dopt_seqPPH.temp = calcEPPHseqdata(y_doptimal, x_doptimal, models, sigmasq, initN, numSeq, N_seq)
  dopt_seqPPH0[ , j] = dopt_seqPPH.temp[1, ]
  dopt_seqPPH1[ , j] = dopt_seqPPH.temp[2, ]
  fact_seqPPH.temp = calcEPPHseqdata(y_factorial, x_factorial, models, sigmasq, initN, numSeq, N_seq)
  fact_seqPPH0[ , j] = fact_seqPPH.temp[1, ]
  fact_seqPPH1[ , j] = fact_seqPPH.temp[2, ]
  smmed_seqPPH.temp = calcEPPHseqdata(smmedsim$y, smmedsim$D, models, sigmasq, initN, numSeq, N_seq)
  smmed_seqPPH0[ , j] = smmed_seqPPH.temp[1, ]
  smmed_seqPPH1[ , j] = smmed_seqPPH.temp[2, ]
}

par(mfrow=c(1,2))

#mean
rand_seqPPH0_mean = apply(rand_seqPPH0, 1, mean)
dopt_seqPPH0_mean = apply(dopt_seqPPH0, 1, mean)
fact_seqPPH0_mean = apply(fact_seqPPH0, 1, mean)
smmed_seqPPH0_mean = apply(smmed_seqPPH0, 1, mean)
rand_seqPPH1_mean = apply(rand_seqPPH1, 1, mean)
dopt_seqPPH1_mean = apply(dopt_seqPPH1, 1, mean)
fact_seqPPH1_mean = apply(fact_seqPPH1, 1, mean)
smmed_seqPPH1_mean = apply(smmed_seqPPH1, 1, mean)
# PPH0
idxlast = numSeq + 1
plot(x = 1:idxlast, y = smmed_seqPPH0_mean, type = "l",
     xlab = "", ylab = "E[P(H0|Y,X)|X]", ylim = c(0, 1))
lines(x = 1:idxlast, rand_seqPPH0_mean, col = 3, lty = 2)
lines(x = 1:idxlast, dopt_seqPPH0_mean, col = 4, lty = 2)
lines(x = 1:idxlast, fact_seqPPH0_mean, col = 6, lty = 2)
points(x = idxlast, y = rand_seqPPH0_mean[idxlast], col = 3, pch = 16)
points(x = idxlast, y = dopt_seqPPH0_mean[idxlast], col = 4, pch = 17)
points(x = idxlast, y = fact_seqPPH0_mean[idxlast], col = 6, pch = 18)
points(x = idxlast, y = smmed_seqPPH0_mean[idxlast])
legend("topright", design_names, lty = c(1, rep(2, length(design_names) - 1)), 
       pch = c(1, 16, 17, 18), col = design_col, bg = "white")
# PPH1
plot(x = 1:idxlast, y = smmed_seqPPH1_mean, type = "l",
     xlab = "", ylab = "E[P(H1|Y,X)|X]", ylim = c(0, 1))
lines(x = 1:idxlast, rand_seqPPH1_mean, col = 3, lty = 2)
lines(x = 1:idxlast, dopt_seqPPH1_mean, col = 4, lty = 2)
lines(x = 1:idxlast, fact_seqPPH1_mean, col = 6, lty = 2)
points(x = idxlast, y = rand_seqPPH1_mean[idxlast], col = 3, pch = 16)
points(x = idxlast, y = dopt_seqPPH1_mean[idxlast], col = 4, pch = 17)
points(x = idxlast, y = fact_seqPPH1_mean[idxlast], col = 6, pch = 18)
points(x = idxlast, y = smmed_seqPPH1_mean[idxlast])



# --- MSE(Bn) --- #

par(mfrow = c(1, length(betaT)))

smmedvs = smmedvsH1
hyp_mu = mu1
hyp_V = V1
hyp_ind = indices1

mseBn_smmed = getClosedMSE(smmedvs$D, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_rand = getClosedMSE(x_random, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_dopt = getClosedMSE(x_doptimal, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_fact = getClosedMSE(x_factorial, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean

for(i in 1:length(betaT)){
  barplot(c(mseBn_smmed[i], mseBn_rand[i], mseBn_dopt[i], mseBn_fact[i]), names.arg = design_names, las = 2)
}



# --- the design --- #

par(mfrow = c(1, 3))

numBreaks = 8
smmed = smmedvsH1
Ntot = dim(smmed$D)[1]

maxcounts = rep(NA, length(mu_full))
for(i in 1:length(mu_full)){
  marginal = i
  h = hist(smmed$D[ (initN + 1):Ntot, marginal], breaks = numBreaks, plot = FALSE)
  maxcounts[i] = max(h$counts)
}

for(i in 1:length(mu_full)){
  marginal = i
  h1 = hist(smmed$D[ (initN + 1):Ntot, marginal], breaks = numBreaks, 
            ylim = range(0, maxcounts), main = "", 
            xlab = paste("marginal ", i, sep = ""))
}

