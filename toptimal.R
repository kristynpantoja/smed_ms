library(rodd)
library(mvtnorm)

sigmasq = 0.5
xmin = -1
xmax = 1
N = 50

betaT = c(-0.2, 0, 0.5) # like ex2.1
typeT = 3

numSims = 500
space_filling = seq(from = xmin, to = xmax, length.out = N)
ysims_space = simulateY(space_filling, N, betaT, sigmasq, numSims = numSims, typeT, seed = 123)

lm_lin_sims = matrix(NA, 2, numSims)
lm_quad_sims = matrix(NA, 3, numSims)
lm_cub_sims = matrix(NA, 4, numSims)
for(i in 1:numSims){
  lm_lin_sims[ , i] = lm(ysims_space[ , i] ~ space_filling)$coefficients
  lm_quad_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2))$coefficients
  lm_cub_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2) + 
                           I(space_filling^3))$coefficients
}
lm_lin = apply(lm_lin_sims, 1, mean)
lm_quad = apply(lm_quad_sims, 1, mean)
lm_cub = apply(lm_cub_sims, 1, mean)





# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)
?rodd

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta <- list(eta0, eta1)

#List of fixed parameters
theta0 <- lm_lin
theta1 <- lm_quad
theta.fix <- list(theta0, theta1)

#Comparison table, holding quadratic (i = 2) model parameters fixed while letting linear model parameters vary
p <- matrix(c(0, 0, 
              1, 0), c(2, 2), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)


#Comparison table, holding linear (i = 1) model parameters fixed while letting quadratic model parameters vary
p2 <- matrix(c(0, 1, 
              0, 0), c(2, 2), byrow = TRUE)
#Design estimation - is close to D-Optimal quadratic design, but with different weights (more in middle), and slightly off-center
res2 <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p2, 
             x.lb = -1, x.rb = 1)








# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta2 <- function(x, theta2) 
  theta2[1] + theta2[2] * x + theta2[3] * x^2 + theta2[4] * x^3

eta <- list(eta0, eta1, eta2)

#List of fixed parameters
theta0 <- lm_lin
theta1 <- lm_quad
theta2 <- lm_cub
theta.fix <- list(theta0, theta1, theta2)

#Comparison table, holding quadratic (i = 2) model parameters fixed while letting cubic model parameters vary
p <- matrix(c(0, 0, 0,
              0, 0, 1,
              0, 0, 0), c(3, 3), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)


#Comparison table, holding cubic (i = 1) model parameters fixed while letting quadratic model parameters vary
p2 <- matrix(c(0, 0, 0,
              0, 0, 0,
              0, 1, 0), c(3, 3), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res2 <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p2, 
             x.lb = -1, x.rb = 1)


















#######################
#######################

# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)
?rodd

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta <- list(eta0, eta1)

#List of fixed parameters
theta0 <- c(1, 1)
theta1 <- c(1, 1, 1)
theta.fix <- list(theta0, theta1)

#Comparison table
p <- matrix(c(0, 0, 
              1, 0), c(2, 2), byrow = TRUE)

#Design estimation - is the D-Optimal design... why?
res <- tpopt(x = c(-1, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)
res <- tpopt(x = c(-1, 0, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)

plot(res)
summary(res)






######
#####
####
###
##
#
##
###
####
#####
######

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.3 # note: needed to change this value to make 

xmin = -1
xmax = 1
N = 100

numSims = 500
space_filling = seq(from = xmin, to = xmax, length.out = N)
ysims_space = simulateY(space_filling, N, betaT, sigmasq, numSims = numSims, typeT, seed = 123)

lm_lin_sims = matrix(NA, 2, numSims)
lm_quad_sims = matrix(NA, 3, numSims)
lm_cub_sims = matrix(NA, 4, numSims)
for(i in 1:numSims){
  lm_lin_sims[ , i] = lm(ysims_space[ , i] ~ space_filling)$coefficients
  lm_quad_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2))$coefficients
  lm_cub_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2) + 
                           I(space_filling^3))$coefficients
}
lm_lin = apply(lm_lin_sims, 1, mean)
lm_quad = apply(lm_quad_sims, 1, mean)
lm_cub = apply(lm_cub_sims, 1, mean)

# install.packages("rodd")
library(rodd)

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta <- list(eta0, eta1)

#List of fixed parameters
theta0 <- lm_lin
theta1 <- lm_quad
theta.fix <- list(theta0, theta1)

#Comparison table, holding quadratic (i = 2) model parameters fixed while letting linear model parameters vary
p <- matrix(c(0, 0, 
              1, 0), c(2, 2), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)
topt1 = c(rep(-1, round(0.25*N)), rep(0, N - 2*round(0.25*N)), rep(1, round(0.25*N)))

#Comparison table, holding linear (i = 1) model parameters fixed while letting quadratic model parameters vary
p2 <- matrix(c(0, 1, 
               0, 0), c(2, 2), byrow = TRUE)
#Design estimation - is close to D-Optimal quadratic design, but with different weights (more in middle), and slightly off-center
res2 <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p2, 
              x.lb = -1, x.rb = 1)
topt2 = c(rep(-1, round(0.2463913*N)), rep(0.1553860, N - round(0.2463913*N)))

# MED design #
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^3
k = 4
p = 2
N = 100
numSeq = 10
N_seq = 10
max_alpha = 2 * p
alpha_seq = c(0, ((2:numSeq) / numSeq) * (max_alpha))
smmed_data = readRDS(paste(home, "/smmed/designs/smmed1.rds", sep = ""))
smmed_data_list = readRDS(paste(home, "/smmed/designs/case1smmeds.rds", sep = ""))
numSeqMMED = length(smmed_data_list)

dopt_linear = c(rep(1,floor(N/2)),rep(-1, N - floor(N/2)))
dopt_quadratic = c(rep(1,floor(N/3)),rep(0,ceiling(N/3)),rep(-1,N - floor(N/3) - ceiling(N/3)))
design_names = c("doptlin", "doptquad", "spacefill", "SMMED")
design_col = c(4, 5, 3, 1)

par(mfrow=c(2,2))
hist(space_filling, breaks = 20, main = "Space-Filling")
hist(smmed_data$D, breaks = 20, main = "Sequential M-MED")
hist(dopt_linear, breaks = 20, main = "Linear D-Optimal")
hist(dopt_quadratic, breaks = 20, main = "Quadratic D-Optimal")

exppostprobs_space = calcExpPostProbH(space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt1 = calcExpPostProbH(dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt2 = calcExpPostProbH(dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_topt1 = calcExpPostProbH(topt1, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_topt2 = calcExpPostProbH(topt2, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)

# calculate posterior probabilities of hypotheses based on this data
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

par(mfrow=c(1,3))

plot(x = 1:numSeq, y = postprobs0, type = "l", main = "P(H0|Y)", xlab = "steps 1:numSeq", ylab = "P(H0|Y)", ylim = range(postprobs0, exppostprobs_space[1], exppostprobs_dopt1[1], exppostprobs_dopt2[1]))
for(k in 1:numSeqMMED){
  lines(x = 1:numSeq, y = postprobs0seq[,k], col=rgb(0, 0, 0, 0.1))
}
lines(x = 1:numSeq, y = postprobs0seq_mean, lty = 2)
abline(h = exppostprobs_space[1], col = 3)
abline(h = exppostprobs_dopt1[1], col = 4)
abline(h = exppostprobs_dopt2[1], col = 5)
abline(h = exppostprobs_topt1[1], col = 6, lty = 2)
abline(h = exppostprobs_topt2[1], col = 7, lty = 2)
legend("topright", c(design_names, "SMMEDsims", "SMMEDsimsmean", "topt1", "topt2"), lty = c(rep(1, length(design_col)), 1, 2, 2, 2), col = c(design_col, rgb(0, 0, 0, 0.1), 1, 6, 7))

plot(x = 1:numSeq, y = postprobs1, type = "l", main = "P(H1|Y)", xlab = "steps 1:numSeq", ylab = "P(H1|Y)", ylim = range(postprobs1, exppostprobs_space[2], exppostprobs_dopt1[2], exppostprobs_dopt2[2]))
for(k in 1:numSeqMMED){
  lines(x = 1:numSeq, y = postprobs1seq[,k], col=rgb(0, 0, 0, 0.1))
}
lines(x = 1:numSeq, y = postprobs1seq_mean, lty = 2)
abline(h = exppostprobs_space[2], col = 3)
abline(h = exppostprobs_dopt1[2], col = 4)
abline(h = exppostprobs_dopt2[2], col = 5)
abline(h = exppostprobs_topt1[2], col = 6, lty = 2)
abline(h = exppostprobs_topt2[2], col = 7, lty = 2)

plot(x = 1:numSeq, y = BF01s, type = "l", main = "BF01", xlab = "steps 1:numSeq", ylab = "BF01", ylim = range(BF01s, exppostprobs_space[3], exppostprobs_dopt1[3], exppostprobs_dopt2[3]))
for(k in 1:numSeqMMED){
  lines(x = 1:numSeq, y = BF01seq[,k], col=rgb(0, 0, 0, 0.1))
}
lines(x = 1:numSeq, y = BF01seq_mean, lty = 2)
abline(h = exppostprobs_space[3], col = 3)
abline(h = exppostprobs_dopt1[3], col = 4)
abline(h = exppostprobs_dopt2[3], col = 5)
abline(h = exppostprobs_topt1[3], col = 6, lty = 2)
abline(h = exppostprobs_topt2[3], col = 7, lty = 2)


