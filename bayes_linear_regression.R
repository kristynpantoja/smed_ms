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
var_mean = 0.05
var_e = 0.1 # same variance

## Running Algorithm

n = 23
numCandidates = 1000
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)





# experimenting with Lattice function

library(OptimalDesign)
?od.m1
candidates = runif(numCandidates, xmin, xmax)
od.m1(candidates, b, A=NULL, w0=NULL, type="exact", kappa=1e-9,
      tab=NULL, graph=NULL, t.max=120)

library(mined)
?Lattice
mined::Lattice(7, 1)

