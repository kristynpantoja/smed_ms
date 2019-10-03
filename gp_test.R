library(MASS)
library(reshape2)
library(ggplot2)
library(transport)
library(mined)
# --- Begin Helper Functions --- #
# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
#C_fn = function(X, Y) outer(X, Y, FUN = function(X, Y) C_fn_elementwise(X, Y, l))
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
#C_fn2 = function(X, Y) outer(X, Y, C_fn_elementwise, l)
# --- End Helper Functions --- #
# --- Begin GP Notes from Example --- #
# GP with Squared Exponential C
# resource for getting started: 
#  https://www.carlboettiger.info/2012/10/17/basic-regression-in-gaussian-processes.html
# GP Parameter Specifications
l = 0.2 # scale parameter for radial covariance function
# observation: the smaller the "l" the more wiggly the plot.
## With Data
# In general we aren’t interested in drawing from the prior,
#  but want to include information from the data as well.
#  Here, we assume no noise.
# data:
obs = data.frame(x = c(-4, -3, -1,  0,  2), y = c(-2,  0,  1,  2, -1))
# new design points
xstar = seq(-5,5,len=50)
# parameter
l = 1
#C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
## draws from null distribution (prior at those points xstar)
# Get covariance matrix from Covariance function
covariance_xstar = C_fn(xstar, xstar, l)
# Generate a number of functions from the GP, mean = 0 for now.
numFns = 2
set.seed(1234)
values <- mvrnorm(numFns, rep(0, length=length(xstar)), covariance_xstar)

# R eshape the data into long (tidy) form, listing x value, y value, and sample number
# Reshape the data into long (tidy) form, listing x value, y value, and sample number
data <- data.frame(x=xstar, t(values))
data <- melt(data, id="x")
head(data)
# Plot
fig1 <- ggplot(data,aes(x=x,y=value)) +
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig1
# Assume y ~ N(0, cov(X, X))
# Then the conditional probability of observing our data is
#  y* | X*, X, y ~ N(K(X*, X) K(X, X)^(-1) y, 
#                    K(X*, X*) - K(X*, X) K(X, X)^(-1) K(X, X*))
# Calculate posterior mean and covariance:
cov_xx_inv <- solve(C_fn(obs$x, obs$x, l))
Ef <- C_fn(xstar, obs$x, l) %*% cov_xx_inv %*% obs$y
Cf <- C_fn(xstar, xstar, l) - C_fn(xstar, obs$x, l) %*% cov_xx_inv %*% C_fn(obs$x, xstar, l)
# Draw random samples from posterior
set.seed(1234)
values <- mvrnorm(3, Ef, Cf)
dat = data.frame(x=xstar, t(values))
dat = melt(dat, id="x")
dat2 = data.frame(xstar, Ef)
fig2 = ggplot(dat,aes(x=dat$x,y=dat$value)) +
  geom_ribbon(data=dat2, 
              aes(x=xstar, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), 
                  ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") +
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=dat2,aes(x=xstar,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2
## Noise in Data
# y = f(x) + epsilon (i.e. some process noise)
# epsilon ~ Normal(0, sigma_eps)
sigma_eps <- 0.25
cov_xx_inv_noise <- solve(C_fn(obs$x, obs$x, l) + sigma_eps^2 * diag(1, length(obs$x)))
E_f_noise <- C_fn(xstar, obs$x, l) %*% cov_xx_inv_noise %*% obs$y
C_f_noise <- C_fn(xstar, xstar, l) - C_fn(xstar, obs$x, l) %*% cov_xx_inv_noise %*% C_fn(obs$x, xstar, l)
# 3 samples from posterior 
set.seed(1234)
values_noise <- mvrnorm(3, E_f_noise, C_f_noise)
# plot
dat_noise <- data.frame(x=xstar, t(values_noise))
dat_noise <- melt(dat_noise, id="x")
dat_noise2 = data.frame(xstar, E_f_noise)
fig3 <- ggplot(dat_noise,aes(x=x,y=value)) +
  geom_ribbon(data=dat_noise2, 
              aes(x=xstar, y=E_f_noise, ymin=(E_f_noise-2*sqrt(diag(C_f_noise))), 
                  ymax=(E_f_noise+2*sqrt(diag(C_f_noise)))),
              fill="grey80") + # Var
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=dat_noise2,aes(x=xstar,y=E_f_noise), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig3
# --- End GP Notes from Example --- #
# --- Begin GP Other Notes / Playing Around --- #
## choose appropriate parameter l and other stuff for points X in [0, 1]
############
## Drawing from GP Prior (i.e. null distribution), for X between 0 and 1
# Try some X
X = seq(from = 0, to = 1, length.out = 51)
# GP Parameter Specifications
l = 0.2 # scale parameter for radial covariance function
# observation: the smaller the "l" the more wiggly the plot.
# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
#C_fn = function(X, Y) outer(X, Y, FUN = function(X, Y) C_fn_elementwise(X, Y, l))
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
#C_fn2 = function(X, Y) outer(X, Y, C_fn_elementwise, l)
# Get covariance matrix from Covariance function
covariance = C_fn(X, X, l)
# Generate a number of functions from the GP, mean = 0 for now.
numFns = 2
set.seed(1234)
values <- mvrnorm(numFns, rep(0, length=length(X)), covariance)
# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data <- data.frame(x=X, t(values))
data <- melt(data, id="x")
head(data)
# Plot
fig1 <- ggplot(data,aes(x=x,y=value)) +
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig1
# this is a plot of draws from the null distribution.
#data1 = subset(data, variable == "X1")
#data2 = subset(data, variable == "X2")
#plot(x = data1$x, data1$value, type = "l")
############
# ...
## GPs Different l 
# Try some X
X = seq(from = 0, to = 1, length.out = 51)
# GP Parameter Specifications
l1 = 0.2 # scale parameter for radial covariance function
l2 = 0.3
# observation: the smaller the "l" the more wiggly the plot.
# Covariance function
# still using radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
# Get covariance matrix from Covariance function
cov1 = C_fn(X, X, l1)
cov2 = C_fn(X, X, l2)
# Generate a fn from each GP
set.seed(1234)
values1 <- mvrnorm(1, rep(0, length=length(X)), cov1)
values2 <- mvrnorm(1, rep(0, length=length(X)), cov2)
# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data1 = data.frame(x = X, values = values1)
data2 = data.frame(x = X, values = values2)
# Plot
ymin = min(c(data1$values, data2$values))
ymax = max(c(data1$values, data2$values))
plot(data1, type = "l", ylim = c(ymin, ymax))
lines(data2)



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################   Gaussian Process Model Selection   #####################
####################################################################################
####################################################################################
####################################################################################
####################################################################################



### Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)

Wasserstein_distance = function(mu1, mu2, var1, var2, type = 1){
  if(type == 1){ # univariate case
    wass = sqrt((mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
  } else{
    if(type > 1){ # multivariate case 
      sqrt_var2 = sqrtm(var2)
      wass = sqrt(crossprod(mu1 - mu2) + sum(diag(var1 + var2 - 2 * sqrtm(sqrt_var2 %*% var1 %*% sqrt_var2))))
    } else{
      stop("invalid type")
    }
  }
  return(as.numeric(wass))
}

## test Wasserstein function

# multivariate case

p = 4
mu1 = rep(7, 4)
var1 = 1 * diag(p)
mu2 = rep(7, 4)
var2 = matrix(0.7, p, p) + diag(rep(1-0.7, p))
N = 500

Wasserstein_distance(mu1, mu2, var1, var2, type = 2) # 1.091925

library(transport)
library(mnormt)
X1 = rmnorm(N, mean = mu1, varcov = var1)
X2 = rmnorm(N, mean = mu2, varcov = var2)
wasserstein(pp(X1), pp(X2), p) # 1.545051
##### is this close enough?? #####

# univariate case
mu1 = 1
mu2 = 20
var1 = 4
var2 = 9
Wasserstein_distance(mu1, mu2, var1, var2, type = 1) # 19.0263
x1 = rnorm(N, mean = mu1, sd = sqrt(var1))
x2 = rnorm(N, mean = mu2, sd = sqrt(var2))
wasserstein1d(x1, x2) # 19.05987
# definitely works here!



## 

q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly(x, var_e, var_mean) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}

f_min_fast = function(candidate_jk, D_k, gamma_k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate_jk, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k * 
    max(sapply(D_k, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean)^gamma_k / abs(x_i - candidate_jk))))
}

avg_w = function(x, D_k, j, N, K0, K1, T){
  Xnotj = D_k[-j]
  numObs = length(D_k[-j])
  if(numObs != (N - 1)) print("length(D_k[-j]) is NOT equal to (N - 1)")
  wassersteins = rep(NA, T)
  Kernel = NULL # different for each simulation t
  which_kernel = round(runif(T))
  y_notj = matrix(rep(NA, numObs * T), nrow = numObs, ncol = T) # each column is a sim. of y_notj
  for(t in 1:T){
    # Generate y
    if(which_kernel[t] == 0){
      Kernel = K0
    } else{
      Kernel = K1
    }
    cov_Xnotj_Xnotj = Kernel(Xnotj, Xnotj)
    cov_xx_inv <- solve(C_fn(obs$x, obs$x, l))
    # Simulate y_notj^(t)
    ysimt = mvrnorm(1, rep(0, length = numObs), cov_Xnotj_Xnotj)
    y_notj[ , t] = ysimt
    # Get predictive means and covariances
    # m = 0:
    cov_Xnotj_Xnotj_inv0 = solve(K0(Xnotj, Xnotj))
    M0x = K0(x, Xnotj) %*% cov_Xnotj_Xnotj_inv0 %*% ysimt
    Sigma0x = K0(x, x) - K0(x, Xnotj) %*% cov_Xnotj_Xnotj_inv0 %*% K0(Xnotj, x)
    # m = 1
    cov_Xnotj_Xnotj_inv1 = solve(K1(Xnotj, Xnotj))
    M1x = K1(x, Xnotj) %*% cov_Xnotj_Xnotj_inv1 %*% ysimt
    Sigma1x = K1(x, x) - K1(x, Xnotj) %*% cov_Xnotj_Xnotj_inv1 %*% K1(Xnotj, x)
    wassersteins[t] = Wasserstein_distance(M0x, M1x, Sigma0x, Sigma1x, type = 2)
  }
  return(mean(wassersteins))
}

SMED_ms_fast = function(k0, k1, N = 11, xmin = 0, xmax = 1, K, p = 1, T = 3){
  # k0 is the covariance function / kernel proposed by H0
  # k1 is the covariance function / kernel proposed by H1
  # each of k0 and k1 take in two arguments: xi, xj
  # N is the number of design points
  
  # K0 and K1 are the covariance matrices that result from the outer product of k0 and k1 resp.
  K0 = function(X, Y) outer(X, Y, FUN = k0)
  K1 = function(X, Y) outer(X, Y, FUN = k1)
  
  # -- Make D1 -- #
  
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = seq(from = xmin, to = xmax, length.out = N)
  D1 = C1
  
  # -- If K = 1, return the space-filling design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose next design -- #
  D = matrix(rep(D1, K), nrow = N, ncol = K)
  gammas = c(1:K) / (K - 1)
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # for j = 1
    # get candidates in neighborhood L_1k = (lower, upper):
    R1k = min(abs(D[-1, k] - D[1, k]))
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0? #######
    L1k_upper = min(D[1, k] + R1k, 1) #IT IS, bc o/w GET NaNs in q evaluation!
    # candidates from space-filling design = tildeD1_kplus1
    tildeD1_kplus1 = seq(from = L1k_lower, to = L1k_upper, length.out = N)
    # update C1_kplus1: save the candidates to be used in future designs
    C[[1]] = c(C[[1]], tildeD1_kplus1)
    # criterion to choose first candidate from candidate set:
    # the point at which f1 and f2 are most different
    q_evals = sapply(C[[1]], FUN = function(x) q(x, mean_beta0, mean_beta1, var_e, var_mean))
    xinitind = which.min(q_evals)
    D[1, k + 1] = C[[1]][xinitind] # x1, first element of new set of design points, D_k+1
    
    # for j = 2:n
    for(j in 2:N){
      if(j == N){
        # get candidates in neighborhood L_jk = (lower, upper)
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = min(D[j, k] + R_jk, 1)
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus1 = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDj_kplus1) # This is now C_j^{k+1}
        # evaluate criterion for each point, choose that with smallest criterion evaluation
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        # get candidates in neighborhood L_jk = (lower, upper)
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDjk) # This is now C_j^{k+1}
        # evaluate criterion for each point, choose that with smallest criterion evaluation
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(x, D[1:(j - 1), k + 1], gammas[k], mean_beta0, mean_beta1, var_e, var_mean))
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return("D" = D, "candidates" = C))
}