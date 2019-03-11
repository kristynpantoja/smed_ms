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


### Helper  Functions ###

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2){
  return(sqrt(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
}

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_e + x^2 * var_mean # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min = function(candidate, D, k, mean_beta0, mean_beta1, var_e, var_mean){
  q(candidate, mean_beta0, mean_beta1, var_e, var_mean)^k * 
    sum(sapply(D, function(x_i) (q(x_i, mean_beta0, mean_beta1, var_e, var_mean) / abs(x_i - candidate))^k))
}

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(mean_beta0, mean_beta1, var_e, var_mean, n = 10, numCandidates = 1000, 
                   k = 4, xmin = 0, xmax = 1){
  #  f0, mean_beta0 : regression line and mean slope of null model
  #  f1, mean_beta1 : regression line and mean slope of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on inputs
  
  # Draw a slope for each model
  beta0 = rnorm(n = 1, mean = mean_beta0, sd = sqrt(var_mean))
  beta1 = rnorm(n = 1, mean = mean_beta1, sd = sqrt(var_mean))
  
  # Create linear model
  f0 = function(x) mean_beta0 * x # null regression model
  f1 = function(x) mean_beta1 * x # alternative regression model
  
  # Calculate posterior means and variance for the mean beta_i given the data y - DON'T NEED THESE.
  # where m = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta_i) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # and s = ((1 / sigma_e^2) + (1 / var_mean^2))^(-1)
  # post_beta0 = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta0) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # post_beta1 = ((1 / sigma_e^2) * y + (1 / var_mean^2) * mean_beta1) / ((1 / sigma_e^2) + (1 / var_mean^2))
  # post_var = ((1 / sigma_e^2) + (1 / var_mean^2))^(-1)
  
  # -- Generate Candidate Points -- #
  candidates = runif(numCandidates, xmin, xmax)
  
  # -- Initialize 1st Design Point in D -- #
  # get the point at which f1 and f2 are most different
  xinitind = which.max(abs(f0(candidates) - f1(candidates)))
  D = candidates[xinitind] # x1, first element of set of design points, D
  # also D is used to initialize greedy algorithm
  
  # candidates that are leftover (options to choose from for next 2:n design points)
  candidates = candidates[-xinitind] # candidate set, for choosing next design point x_{n+1}
  
  # Plot density and (highest lik) points
  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
  curve(f1, col = 1, add = TRUE)
  text(D, f0(D), 1, col=4)
  text(D, f1(D), 1, col=4)
  points(D, 0, col=2)
  
  # Sequentially pick rest
  for(i in 2:n){
    # Find f_opt: minimum of f_min
    
    #f_opt = which.min(f_min(candidates, D, k, mean_beta0, mean_beta1, var_e, var_mean))
    f_min_candidates = sapply(candidates, function(x) f_min(x, D, k, mean_beta0, mean_beta1, var_e, var_mean))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    # Update set of design points (D) and plot new point
    D = c(D,xnew) # add the new point to the set
    #candidates = candidates[-f_opt]
    text(xnew, f0(xnew), i, col = 4)
    text(xnew, f1(xnew), i, col = 4)
    points(xnew, 0, col = 2)
  }
  return(D)
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.01
var_e = 1 # same variance

n = 50
numCandidates = 1000
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_linear_regression.png')
#dev.off()




# pictures

n = 50
numCandidates = 1000
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k4.png')
#dev.off()



n = 50
numCandidates = 1000
k = 100
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k100.png')
#dev.off()



n = 50
numCandidates = 1000
k = 125
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k125.png')
#dev.off()



n = 50
numCandidates = 1000
k = 150
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k150.png')
#dev.off()



n = 50
numCandidates = 1000
k = 200
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k200.png')
#dev.off()



n = 50
numCandidates = 1000
k = 500
xmin = 0
xmax = 1

X_test = SMED_ms(mean_beta0, mean_beta1, var_e, var_mean, n, numCandidates, k, xmin, xmax)

#dev.copy(png,'bayes_lr_k500.png')
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

