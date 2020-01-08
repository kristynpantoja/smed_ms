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

# --- sources to designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/postmean_mse_mc.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
source(paste(functions_home, "/plot_EPH1.R", sep = ""))
source(paste(functions_home, "/plot_mse.R", sep = ""))

# choose example: 1, 2, 2.1, 3, 4, 4.1, 5
example = 4.1

# shared parameters/settings:
sigmasq01 = 0.01
sigmasq = 1
# function settings (including and based on prior settings above)
numCandidates = 10^3
k = 4
S = 10
xmin = -1
xmax = 1
p = 2
N = 100

seed = 123
numSims = 10
N = 100

if(example == 1){
  # --- example 1 --- #
  
  # MED design #
  # Priors
  mu0 = c(0,1)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(0,1,0)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex1.RData", sep = ""))
  designs = designs1d_ex1
  
  typeT = 3
  betaT = c(-5,1,5) #mu1+10*rep(sigmasq01,length(mu1))
}
if(example == 2){
  # --- example 2 --- #
  
  # MED design #
  # Priors
  mu0 = c(0, 0)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(-0.2, 0, 0.4)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex2.RData", sep = ""))
  designs = designs1d_ex2
  
  typeT = 3
  betaT = mu1
}
if(example == 2.1){
  # --- example 2 --- #
  
  # MED design #
  # Priors
  mu0 = c(0, 0)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(-0.2, 0, 0.5)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex2pt1.RData", sep = ""))
  designs = designs1d_ex2.1
  
  typeT = 3
  betaT = mu1
}
if(example == 3){
  # --- example 3 --- #
  
  # MED design #
  # Priors
  mu0 = c(0, 0)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(1, 0, -1)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex3.RData", sep = ""))
  designs = designs1d_ex3
  
  typeT = 3
  betaT = mu1
}
if(example == 4){
  # --- example 4 --- #
  
  # MED design #
  # Priors
  mu0 = c(1, 1)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(1, 0, -1)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex4.RData", sep = ""))
  designs = designs1d_ex4
  
  typeT = 3
  betaT = mu1
}
if(example == 4.1){
  # --- example 4.1 --- #
  
  # MED design #
  # Priors
  mu0 = c(0.2, 0.2)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(0.2, 0, -0.2)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
  type01 = c(2, 3)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex4pt1.RData", sep = ""))
  designs = designs1d_ex4.1
  
  typeT = 3
  betaT = mu1
}
if(example == 5){
  # --- example 4 --- #
  
  # MED design #
  # Priors
  mu0 = c(0, 0)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(0, 0.2)
  V1 = diag(rep(sigmasq01,length(mu1)))
  # function settings
  f0 = function(x) mu0[1] + mu0[2] * x
  f1 = function(x) mu1[1] + mu1[2] * x
  type01 = c(2, 2)
  
  # load designs
  load(paste(home, "/designs/designs1d_ex5.RData", sep = ""))
  designs = designs1d_ex5
  
  typeT = 2
  betaT = mu1
}

# # see the models (H0, H1, true)
# par(mfrow=c(2,2))
# mu0.test =mu0
# mu1.test = mu1
# betaT.test = betaT
# par(mfrow = c(2,2))
# xs = seq(-1, 1, length.out = 51)
# modelH0 = function(x){
#   desX = matrix(c(1, x))
#   f = as.vector(t(desX) %*% mu0.test)
#   return(f)
# }
# ysH0 = sapply(xs, FUN = modelH0)
# plot(xs, ysH0, xlim = c(-1, 1), ylim = c(-1, 2), main = "H0 model", type ="l")
# modelH1 = function(x){
#   desX = matrix(c(1, x, x^2))
#   f = as.vector(t(desX) %*% mu1.test)
#   return(f)
# }
# ysH1 = sapply(xs, FUN = modelH1)
# plot(xs, ysH1, xlim = c(-1, 1), ylim = c(-1, 2), main = "H1 model", type ="l")
# model_true = function(x){
#   desX = matrix(c(1, x, x^2))
#   f = as.vector(t(desX) %*% betaT.test)
#   return(f)
# }
# ys_true = sapply(xs, FUN = model_true)
# plot(xs, ys_true, xlim = c(-1, 1), ylim = c(-1, 2), main = "true model", type ="l")
# lines(xs, ysH0, col = 2); lines(xs, ysH1, col = 4)
# plot(xs, ys_true, main = "true model, full view", type ="l")

par(mfrow=c(1,2))
hist(designs$design_MED_alpha2p, main = "MED,alpha=2p")
hist(designs$design_MED_alpha1, main = "MED,alpha=1")

# see evaluation plots
par(mfrow=c(1,1))
sigmasq01_seq = seq(from = 0, to = 1, length.out = 51)
plot_EPH1(designs, N, betaT, mu0, mu1, sigmasq01_seq, sigmasq, numSims, typeT, type01, seed = seed)
legend("bottomleft", names(designs), col = 1:length(designs), lty = 1)


par(mfrow=c(2, 2))
which.plots = 1:length(mu1)
sigmasq01_seq2 = seq(from = 0, to = 1, length.out = 51)
for(k in 1:length(which.plots)){
  plot_MSE(designs, N, mu1, mu1, sigmasq01_seq, sigmasq, typeT, which.plots = which.plots[k], cols = NULL)
  #legend("topleft", names(designs1d), col = 1:length(designs), lty = 1)
  title(paste("BT","=(", paste(round(mu1,3), collapse = ","), ")", ", MSE(B", k - 1, ")", sep = ""))
}




