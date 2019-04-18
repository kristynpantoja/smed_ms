############################################################
## Based on One-At-a-Time Algorithm, Joseph et. al. 2015 
############################################################

### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

# Null Model: y_i = x_i beta0 + epsilon_i
# vs. 
# Alternative Model: y_i = x_i beta1 + epsilon_i

# where, in both models,
#  epsilon_i ~ N(0, sigma_e^2)
#  y_i ~ N(x_i beta, sigma_e^2), beta = {beta0, beta1}

# but now we put a prior on beta:
# beta_i ~ N(mean_beta_i, var_mean)

# and now we want the distances between the means' posterior distributions
# beta_i | y, x ~ N(m, s)
# where m = ((1 / sigma_e^2) * y + (1 / sigma_mean^2) * mean_beta_i) / ((1 / sigma_e^2) + (1 / sigma_mean^2))
# and s = ((1 / sigma_e^2) + (1 / sigma_mean^2))^(-1)

# note that this can be extended to the case of multiple observed y at a given x, in which case:
# beta_i | y, x ~ N(m, s)
# where m = ((n / sigma_e^2) * y + (1 / sigma_mean^2) * mean_beta_i) / ((n / sigma_e^2) + (1 / sigma_mean^2))
# and s = ((n / sigma_e^2) + (1 / sigma_mean^2))^(-1)

# The marginal of y | x, Hi is no longer N(beta_i * x, sigma_e^2), since now we beta_i has a (prior) distribution
# hence, using iterated expectation and variance,
# y | x, Hi ~ N(mean_beta_i * x, sigma_e^2 + x^2 * sigma_mean^2)

library(transport)
source("smed_ms_functions.R")




### Testing

## Pick Parameters

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.001
var_e = 0.01 # same variance

## Running Algorithm

n = 67
numCandidates = 10^5 # suggested 10^5
k = 4
xmin = 0
xmax = 1

# TIME IT!
ptm <- proc.time()
X_test_1atatime = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)
# when numCandidates = 1500, elapsed = 16.039 (w/o printing plot)
# when numCandidates = 2500, elapsed = 41.744
# when numCandidates = 10^5, elapsed = 1042.209 (w/o printing plot), 1079.530 the second time
proc.time() - ptm

# HISTOGRAM.
hist(X_test_1atatime, freq = T)
mean(X_test_1atatime) # 0.6882634
sd(X_test_1atatime) # 0.2182556

#dev.copy(png, 'oneatatime_seq_pattern2_N67_nC100000_hist.png')
#dev.off()

### --- Distances between x's and y's N = 67 --- ###

Xsort = sort(X_test_1atatime)
Y0 = f0(Xsort)
Y1 = f1(Xsort)
diffX = Xsort[-1] - Xsort[-length(Xsort)]
diffY = Y0 - Y1

sum(diffX) # 0.8859489
sum(diffY) # 23.05683

hist(diffX) # right-skewed, like dist of X
hist(diffY) # left-skewed, but less severe drop off

mean(diffX) # 0.01342347
mean(diffY) # 0.3441317

sd(diffX) # 0.01187598
sd(diffY) # 0.1091278


# PICTURES.
#dev.copy(png, 'oneatatime_N67.png')
#dev.off()


# PICTURES.

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 3))
curve(f1, col = 1, add = TRUE)
y = rep(NA, n)
for(i in 1:n){
  y[i] = i * 0.04
  text(X_test_1atatime[i], y[i], i, col=4)
}
points(X_test_1atatime, rep(0, n), col = 2, pch = "*")
lines(X_test_1atatime, y, col = 3)

#dev.copy(png, 'oneatatime_seq_pattern2_N67_nC100000.png')
#dev.off()


# experimenting with Lattice function

library(OptimalDesign)
?od.m1
candidates = runif(numCandidates, xmin, xmax)
od.m1(candidates, b, A=NULL, w0=NULL, type="exact", kappa=1e-9,
      tab=NULL, graph=NULL, t.max=120)

library(mined)
?Lattice
mined::Lattice(7, 1)















n = 100
numCandidates = 1000
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png, 'oneatatime_N67.png')
#dev.off()

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 4))
curve(f1, col = 1, add = TRUE)

for(i in 1:n){
  text(X_test[i], f1(X_test[i]) + i * 0.035, i, col=4)
}

points(X_k, rep(0, N), col = 2)
