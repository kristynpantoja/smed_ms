# libraries
library(Matrix)
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# source files for evaluations

home = "/home/kristyn/Documents/smed_ms"

# --- Sources for S/MMED --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MMED_nodata.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/add_MMED.R", sep = ""))
source(paste(functions_home, "/SMMED.R", sep = ""))

# --- Sources to evaluate designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/posterior_mean_mse.R", sep = ""))
source(paste(functions_home, "/plot_utils.R", sep = ""))
source(paste(functions_home, "/predictive_yhat_mse.R", sep = ""))

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2

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
p = 1
N = 100

numSeq = 10
N_seq = 10
alpha_seq = 1



# --- other designs --- #

space_filling =  seq(from = xmin, to = xmax, length.out = N)
# doptimal for linear model
dopt_linear = c(rep(1,floor(N/2)),rep(-1, N - floor(N/2)))
# doptimal for quadratic model
dopt_quadratic = c(rep(1,floor(N/3)),rep(0,ceiling(N/3)),rep(-1,N - floor(N/3) - ceiling(N/3)))

#############################
# --- case 1: quadratic --- #
#############################

smmed_data = readRDS(paste(home, "/run_designs/smmed/designs/smmed1.rds", sep = ""))
smmed_data_list = readRDS(paste(home, "/run_designs/smmed/designs/case1smmeds.rds", sep = ""))
numSeqMMED = length(smmed_data_list)

par(mfrow = c(1,2))
hist(smmed_data$D, breaks = 20, main = "", xlab = "x")
plot(x = smmed_data$D, y = smmed_data$y, xlab = "x", ylab = "", col = 1)
legend("bottomleft", legend = c("data", "true model"), lty = c(NA, 1), pch = c(1, NA))
curve(fT, add = TRUE)

# png("image.png", width = 800, height = 600)
# plot(...)
# dev.off()

# --- Fits --- #

designs = list(dopt_linear, dopt_quadratic, space_filling, smmed_data$D)
design_names = c("doptlin", "doptquad", "grid", "smmed")
design_col = c(4, 5, 3, 1)

par(mfrow=c(2,2))

col_postmeans = c(3, 5, 4)
plot(x = 1:numSeq, y = smmed_data$postmean1[1, ], type = "l", 
     ylim = c(-0.5, 0.5), #main = "Posterior Mean, Quadratic", 
     xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
lines(x = 1:numSeq, y = smmed_data$postmean1[2, ], type = "l", col = col_postmeans[2])
lines(x = 1:numSeq, y = smmed_data$postmean1[3, ], type = "l", col = col_postmeans[3])
abline(h = betaT[1], lty = 2, col = col_postmeans[1])
abline(h = betaT[2], lty = 2, col = col_postmeans[2])
abline(h = betaT[3], lty = 2, col = col_postmeans[3])
legend("bottomleft", c("Bn0", "Bn1", "Bn2"), lty = c(1, 1, 1), col = col_postmeans, bg = "white")

col_postmeans = c(8, 6)
plot(x = 1:numSeq, y = smmed_data$postmean0[1, ], type = "l", 
     ylim = c(-0.5, 0.5), #main = "Posterior Mean, Linear", 
     xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
lines(x = 1:numSeq, y = smmed_data$postmean0[2, ], type = "l", col = col_postmeans[2])
legend("topleft", c("Bn0", "Bn1"), lty = c(1, 1), col = col_postmeans, bg = "white")

fEst = function(x) smmed_data$postmean1[1, N_seq] + 
  smmed_data$postmean1[2, N_seq] * x + smmed_data$postmean1[3, N_seq] * x^2
curve(fT, xlim = c(-1, 1))
curve(fEst, col = 2, add = T)
legend("topright", c("true model", "estimated model"), lty = c(1, 1), col = c(1, 2))

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
fEst = function(x) smmed_data$postmean0[1, N_seq] + 
  smmed_data$postmean0[2, N_seq] * x
curve(fT, xlim = c(-1, 1))
curve(fEst, col = 2, add = T)

# wasserstein distance
par(mfrow = c(1, 1))

testx = seq(from = xmin, to = xmax, length.out = numCandidates)
is_already_in_D = testx %in% smmed_data$D
testx = testx[!is_already_in_D]
wass_testx = sapply(testx, function(x) Wasserstein_distance_postpred(x, smmed_data$postmean0[,10], smmed_data$postmean1[,10], 
                                                                     diag(smmed_data$postvar0[,10]), diag(smmed_data$postvar1[,10]), sigmasq, type01))
plot(x = testx, y = (wass_testx), type = "l", ylim = range(-0.5, 0.5), 
     ylab = "wasserstein(x)", xlab = "", main = "")
f0data = function(x) smmed_data$postmean0[,10][1] + smmed_data$postmean0[,10][2] * x
f1data = function(x) smmed_data$postmean1[,10][1] + smmed_data$postmean1[,10][2] * x + smmed_data$postmean1[,10][3] * x^2
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
curve(f0data, col = 2, lwd = 5, add = T)
curve(f1data, add = T, col = 5, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)

# left shade
xshade = seq(from = -1, to = -0.6, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

# middle shade
xshade = seq(from = -0.6, to = 0.6, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

# right shade
xshade = seq(from = 0.6, to = 1, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

legend("bottomleft", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))


# --- EPPH --- #

exppostprobs_space = calcExpPostProbH(space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt1 = calcExpPostProbH(dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt2 = calcExpPostProbH(dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)

changing_postprobs = list()
postprobs0 = rep(NA, numSeq)
postprobs1 = rep(NA, numSeq)
BF01s = rep(NA, numSeq)
for(i in 1:numSeq){
  changing_postprobs[[i]] = calcExpPostProbH_data(smmed_data$y[1:(N_seq * i)], 
                                                  smmed_data$D[1:(N_seq * i)], 
                                                  N = N_seq * i, mu0, V0, mu1, V1, sigmasq, type01)
  postprobs0[i] = changing_postprobs[[i]][1]
  postprobs1[i] = changing_postprobs[[i]][2]
  BF01s[i] = changing_postprobs[[i]][3]
}

# calculate posterior probabilities of hypotheses based on this data
postprobs0seq = matrix(NA, numSeq, numSeqMMED)
postprobs1seq = matrix(NA, numSeq, numSeqMMED)
BF01seq = matrix(NA, numSeq, numSeqMMED)
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcExpPostProbH_data(smmed_data_k$y[1:(N_seq * i)], 
                                               smmed_data_k$D[1:(N_seq * i)], 
                                               N = N_seq * i, mu0, V0, mu1, V1, sigmasq, type01)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    BF01seq[i, k] = changing_postprobs[3]
  }
}

postprobs0seq_mean = apply(postprobs0seq, 1, mean)
postprobs1seq_mean = apply(postprobs1seq, 1, mean)
BF01seq_mean = apply(BF01seq, 1, mean)

par(mfrow=c(1,2))

plot(x = 1:numSeq, y = postprobs0seq_mean, type = "l", main = "", 
     xlab = "steps 1:numSeq", ylab = "", ylim = c(0, 1))
# for(k in 1:numSeqMMED){
#   lines(x = 1:numSeq, y = postprobs0seq[,k], col=rgb(0, 0, 0, 0.1))
# }
# lines(x = 1:numSeq, y = postprobs0seq_mean, lty = 2)
abline(h = exppostprobs_space[1], col = 3)
abline(h = exppostprobs_dopt1[1], col = 4)
abline(h = exppostprobs_dopt2[1], col = 5)
legend("topleft", c(design_names), lty = c(rep(1, length(design_col))), col = c(design_col), bg = "white")

plot(x = 1:numSeq, y = postprobs1seq_mean, type = "l", main = "", 
     xlab = "steps 1:numSeq", ylab = "", ylim = c(0, 1))
# for(k in 1:numSeqMMED){
#   lines(x = 1:numSeq, y = postprobs1seq[,k], col=rgb(0, 0, 0, 0.1))
# }
# lines(x = 1:numSeq, y = postprobs1seq_mean, lty = 2)
abline(h = exppostprobs_space[2], col = 3)
abline(h = exppostprobs_dopt1[2], col = 4)
abline(h = exppostprobs_dopt2[2], col = 5)



# --- MSE(Bn) --- #

postvar_doptlin = getClosedMSE(dopt_linear, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_doptquad = getClosedMSE(dopt_quadratic, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_space = getClosedMSE(space_filling, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_smmed = getClosedMSE(smmed_data$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
par(mfrow = c(1,3))

barplot(c(postvar_doptlin[1], postvar_doptquad[1], postvar_space[1], postvar_smmed[1]), 
        names.arg = design_names, main = "B0")
barplot(c(postvar_doptlin[2], postvar_doptquad[2], postvar_space[2], postvar_smmed[2]), 
        names.arg = design_names, main = "B1")
barplot(c(postvar_doptlin[3], postvar_doptquad[3], postvar_space[3], postvar_smmed[3]), 
        names.arg = design_names, main = "B2")



# --- MSE(y-hat) --- #

x_seq = seq(from = -1, to = 1, length.out = 1e3)
yhatmse_space = getClosedMSEyhat_seq(x_seq, space_filling, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptquad = getClosedMSEyhat_seq(x_seq, dopt_quadratic, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_smmed = getClosedMSEyhat_seq(x_seq, smmed_data$D, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptlin = getClosedMSEyhat_seq(x_seq, dopt_linear, N, betaT, typeT, mu1, V1, sigmasq, type01[2])

par(mfrow = c(1,1))
ylimarg = range(0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_smmed$MSEyhat)
plot(x_seq, yhatmse_space$MSEyhat, type = "l", col = 3, ylim = ylimarg, 
     xlab = "x in [-1, 1]", ylab = "MSE(y-hat)")
lines(x_seq, yhatmse_doptquad$MSEyhat, col = 5)
lines(x_seq, yhatmse_smmed$MSEyhat, col = 1)
lines(x_seq, yhatmse_doptlin$MSEyhat, col = 4)
legend("top", design_names, col = design_col, lty = 1)


#############################
# --- case 2: cubic ------- #
#############################

# MED design #
sigmasq01 = 0.25
sigmasq = 0.1

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
betaT = c(0, -0.75, 0, 1)
typeT = 4
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
xmin = -1
xmax = 1
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
type01 = c(2, 3)
numCandidates = 10^3
k = 4
p = 1
N = 100

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3

numSeq = 10 
N_seq = 10
alpha_seq = 1
smmed_data2 = readRDS(paste(home, "/run_designs/smmed/designs/smmed2.rds", sep = ""))
smmed_data2_list = readRDS(paste(home, "/run_designs/smmed/designs/case2smmeds.rds", sep = ""))
numSeqMMED = length(smmed_data2_list)

par(mfrow = c(1,2))
hist(smmed_data2$D, breaks = 20, main = "", xlab = "x")
plot(x = smmed_data2$D, y = smmed_data2$y, xlab = "x", ylab = "")
legend("bottomleft", legend = c("data", "true model"), lty = c(NA, 1), pch = c(1, NA))
curve(fT, add = TRUE)



# --- fits --- #

par(mfrow=c(2,3))

col_postmeans = c(1, 2, 3, 4)
smmed_data2.postmean2 = matrix(NA, length(betaT), numSeq)
for(k in 1:numSeq){
  smmed_data2.postmean2[ , k] = as.vector(postmean(smmed_data2$y[1:(N_seq * k)], 
                                                   smmed_data2$D[1:(N_seq * k)], (N_seq * k), 
                                                   c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4))
}
plot(x = 1:numSeq, y = smmed_data2.postmean2[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
lines(x = 1:numSeq, y = smmed_data2.postmean2[2, ], type = "l", col = col_postmeans[2])
lines(x = 1:numSeq, y = smmed_data2.postmean2[3, ], type = "l", col = col_postmeans[3])
lines(x = 1:numSeq, y = smmed_data2.postmean2[4, ], type = "l", col = col_postmeans[4])
abline(h = betaT[1], lty = 1, col = col_postmeans[1]) # even though i said lty = 2, it gets covered.
abline(h = betaT[2], lty = 2, col = col_postmeans[2])
abline(h = betaT[3], lty = 2, col = col_postmeans[3])
abline(h = betaT[4], lty = 2, col = col_postmeans[4])
legend("bottomleft", c("Bn0", "Bn1", "Bn2", "Bn3"), lty = 1, 
       col = col_postmeans, bg = "white")

col_postmeans = c(1, 2, 3)
plot(x = 1:numSeq, y = smmed_data2$postmean1[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
lines(x = 1:numSeq, y = smmed_data2$postmean1[2, ], type = "l", col = col_postmeans[2])
lines(x = 1:numSeq, y = smmed_data2$postmean1[3, ], type = "l", col = col_postmeans[3])
legend("bottomleft", c("Bn0", "Bn1", "Bn2"), lty = c(1, 1, 1), col = col_postmeans, bg = "white")

col_postmeans = c(1, 2)
plot(x = 1:numSeq, y = smmed_data2$postmean0[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# would plot error bars from posterior variances, but they're soooo small... why?
lines(x = 1:numSeq, y = smmed_data2$postmean0[2, ], type = "l", col = col_postmeans[2])
legend("bottomleft", c("Bn0", "Bn1"), lty = c(1, 1), col = col_postmeans, bg = "white")

# smmed_data2.postmean2.1 = postmean(smmed_data2$y, smmed_data2$D, N, 
#                                        c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
# fEst = function(x) smmed_data2.postmean2.1[1] + smmed_data2.postmean2.1[2] * x + smmed_data2.postmean2.1[3] * x^2 + smmed_data2.postmean2.1[4] * x^3
fEst = function(x) smmed_data2.postmean2[1, N_seq] + smmed_data2.postmean2[2, N_seq] * x + smmed_data2.postmean2[3, N_seq] * x^2 + smmed_data2.postmean2[4, N_seq] * x^3
curve(fT, xlim = c(-1, 1))
curve(fEst, col = 2, add = T)
legend("bottomleft", c("true model", "estimated model"), lty = c(1, 1), col = c(1, 2))

fEst = function(x) smmed_data2$postmean1[1, N_seq] + smmed_data2$postmean1[2, N_seq] * x + smmed_data2$postmean1[3, N_seq] * x^2
curve(fT, xlim = c(-1, 1))
curve(fEst, col = 2, add = T)

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
fEst = function(x) smmed_data2$postmean0[1, N_seq] + smmed_data2$postmean0[2, N_seq] * x
curve(fT, xlim = c(-1, 1))
curve(fEst, col = 2, add = T)



# --- high density areas (wasserstein distance) --- #

par(mfrow = c(1, 1))
testx = seq(from = xmin, to = xmax, length.out = numCandidates)
is_already_in_D = testx %in% smmed_data$D
testx = testx[!is_already_in_D]
wass_testx = sapply(testx, function(x) Wasserstein_distance_postpred(x, smmed_data$postmean0[,10], smmed_data$postmean1[,10], 
                                                                     diag(smmed_data$postvar0[,10]), diag(smmed_data$postvar1[,10]), sigmasq, type01))
plot(x = testx, y = (wass_testx), type = "l", ylim = range(-1, 1), 
     main = "wasserstein(x)")
f0data = function(x) smmed_data$postmean0[,10][1] + smmed_data$postmean0[,10][2] * x
f1data = function(x) smmed_data$postmean1[,10][1] + smmed_data$postmean1[,10][2] * x + smmed_data$postmean1[,10][3] * x^2
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
curve(f0data, col = 2, lwd = 5, add = T)
curve(f1data, add = T, col = 5, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)

# left shade
xshade = seq(from = -1, to = -0.6, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

# middle shade
xshade = seq(from = -0.6, to = 0.6, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

# right shade
xshade = seq(from = 0.6, to = 1, length.out = 100)
yshade1 = sapply(xshade, FUN = f0data)
yshade2 = sapply(xshade, FUN = f1data)
polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)

legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))



# --- EPPH --- #

designs = list(dopt_linear, dopt_quadratic, space_filling, smmed_data2$D)
design_names = c("doptlin", "doptquad", "spacefill", "SMMED")
design_col = c(4, 5, 3, 1)

# what are expected posterior probabilities of hypotheses for space-filling design? d-optimal designs? 
models = list("H0" = list(mu0, V0, 2),
              "H1" = list(mu1, V1, 3),
              "H2" = list(betaT, diag(rep(sigmasq01, 4)), 4))
numSims = 100
exppostprobs_space = calcEPPH(space_filling, N, betaT, typeT, models, sigmasq, numSims, seed = 123)
exppostprobs_dopt1 = calcEPPH(dopt_linear,  N, betaT, typeT, models, sigmasq, numSims, seed = 123)
exppostprobs_dopt2 = calcEPPH(dopt_quadratic, N, betaT, typeT, models, sigmasq, numSims, seed = 123)

# calculate posterior probabilities of hypotheses based on this data
changing_postprobs = list() 
postprobs = matrix(NA, length(models), numSeq)
for(i in 1:numSeq){
  changing_postprobs[[i]] = calcEPPHdata(smmed_data2$y[1:(N_seq * i)],
                                                          smmed_data2$D[1:(N_seq * i)], 
                                                          N = N_seq * i, models, sigmasq)
  postprobs[ , i] = changing_postprobs[[i]]
}

postprobs0seq = matrix(NA, numSeq, numSeqMMED)
postprobs1seq = matrix(NA, numSeq, numSeqMMED)
postprobs2seq = matrix(NA, numSeq, numSeqMMED)
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data2_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcEPPHdata(smmed_data_k$y[1:(N_seq * i)], 
                                                       smmed_data_k$D[1:(N_seq * i)], 
                                                       N = N_seq * i, models, sigmasq)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    postprobs2seq[i, k] = changing_postprobs[3]
  }
}

postprobs0seq_mean = apply(postprobs0seq, 1, mean)
postprobs1seq_mean = apply(postprobs1seq, 1, mean)
postprobs2seq_mean = apply(postprobs2seq, 1, mean)


numModels = length(models)
postprobs_seq = array(NA, dim = c(numModels, numSeq, numSeqMMED))
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data2_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcEPPHdata(smmed_data_k$y[1:(N_seq * i)], 
                                                       smmed_data_k$D[1:(N_seq * i)], 
                                                       N = N_seq * i, models, sigmasq)
    postprobs_seq[ , i, k] = changing_postprobs
  }
}
postprobs_means = apply(postprobs_seq, c(1,2), mean)

par(mfrow=c(1,length(models)))
for(k in 1:length(models)){
  plot(x = 1:numSeq, y = postprobs_means[ k, ], type = "l", 
       xlab = "steps 1:numSeq", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
       ylim = range(postprobs[ k, ], exppostprobs_space[k], exppostprobs_dopt1[k], exppostprobs_dopt2[k]))
  # for(m in 1:numSeqMMED){
  #   lines(x = 1:numSeq, y = postprobs_seq[k, , m], col=rgb(0, 0, 0, 0.1))
  # }
  # lines(x = 1:numSeq, y = postprobs_means[k, ], col = 1, lty = 2)
  abline(h = exppostprobs_space[k], col = 3)
  abline(h = exppostprobs_dopt1[k], col = 4)
  abline(h = exppostprobs_dopt2[k], col = 5)
  if(k == 1) legend("topright", c(design_names), lty = c(rep(1, length(design_col))), col = c(design_col), bg = "white")
}

# --- MSE(y-hat) --- #

x_seq = seq(from = -1, to = 1, length.out = 1e3)
yhatmse_space = getClosedMSEyhat_seq(x_seq, space_filling, N, betaT, typeT, 
                                     c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_doptquad = getClosedMSEyhat_seq(x_seq, dopt_quadratic, N, betaT, typeT, 
                                        c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_smmed = getClosedMSEyhat_seq(x_seq, smmed_data2$D, N, betaT, typeT, 
                                     c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_doptlin = getClosedMSEyhat_seq(x_seq, dopt_linear, N, betaT, typeT, 
                                       c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
par(mfrow = c(1,1))
ylimarg = range(0, yhatmse_space$MSEyhat, yhatmse_smmed$MSEyhat)
plot(x_seq, yhatmse_space$MSEyhat, type = "l", col = 3, ylim = ylimarg, 
     ylab = "", xlab = "")
lines(x_seq, yhatmse_doptquad$MSEyhat, col = 5)
lines(x_seq, yhatmse_smmed$MSEyhat, col = 1)
lines(x_seq, yhatmse_doptlin$MSEyhat, col = 4)
legend("top", design_names, col = design_col, lty = 1, bg = "white")


