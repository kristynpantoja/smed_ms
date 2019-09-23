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

# choose example: 1, 2

example = 2
if(example == 1){
  # --- example 1 --- #
  
  # MED design #
  # Priors
  sigmasq01 = 0.01
  mu0 = c(0,1)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(0,1,0)
  V1 = diag(rep(sigmasq01,length(mu1)))
  sigmasq = 1
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
  N = 100
  # load designs
  load(paste(home, "/designs/designs1d_ex1.RData", sep = ""))
  designs = designs1d_ex1
  
  # see how values differ with different sigmasq01
  seed = 123
  numSims = 50
  N = 100
  betaT = c(-5,1,5) #mu1+10*rep(sigmasq01,length(mu1))
  type01 = c(2, 3)
  typeT = 3
  sigmasq01_seq = seq(from = 0, to = 0.1, length.out = 51)
}
if(example == 2){
  # --- example 2 --- #
  
  # MED design #
  # Priors, from dj_simulations.R
  sigmasq01 = 0.01
  mu0 = c(0, 0)
  V0 = diag(rep(sigmasq01,length(mu0)))
  mu1 = c(-0.2, 0, 0.4)
  V1 = diag(rep(sigmasq01,length(mu1)))
  sigmasq = 1
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
  N = 100
  # load designs
  load(paste(home, "/designs/designs1d_ex2.RData", sep = ""))
  designs = designs1d_ex2
  
  # see how values differ with different sigmasq01
  seed = 123
  numSims = 100 # 50
  N = 100
  betaT = mu1
  type01 = c(2, 3)
  typeT = 3
  sigmasq01_seq = seq(from = 0, to = 0.1, length.out = 51)
}
plot_EPH1(designs, N, betaT, mu0, mu1, sigmasq01_seq, sigmasq, numSims, typeT, type01, seed = seed)

