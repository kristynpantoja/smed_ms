### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

# Null Model: y_i = x_i beta0 + epsilon_i
# vs. 
# Alternative Model: y_i = x_i beta1 + epsilon_i

# where, in both models,
#  epsilon_i ~ N(0, sigma_e^2)
#  y_i ~ N(x_i beta, sigma_e^2), beta = {beta0, beta1}




library(transport)

### Helper  Functions ###

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2){
  return(sqrt(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
}

# charge function at design point x
q = function(x, beta0, beta1, var_e){
  mu1 = beta0 * x # mean of marginal dist of y | H0
  mu2 = beta1 * x # mean of marginal dist of y | H1
  var = var_e
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min = function(candidate, D, k, beta0, beta1, var_e){
  q(candidate, beta0, beta1, var_e)^k * 
    sum(sapply(D, function(x_i) (q(x_i, beta0, beta1, var_e) / abs(x_i - candidate))^k))
}




### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(f0, beta0, f1, beta1, var_e, n = 10, numCandidates = 1000, k = 4, xmin = 0, xmax = 1, 
                   randomCandidates = T, initialization = 0
                   ){
  #  f0, beta0 : regression line and slope of null model
  #  f1, beta1 : regression line and slope of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on inputs
  # randomCandidates : generate candidate points using uniform distribution by default.
  #                    if FALSE, choose equally-spaced out sequence of candidates
  # initialization : method of initializing first design point of D
  #                  0 = candidate that outputs the largest absolute difference between f0 and f1
  #                  1 = candidate that has the largest ratio of likelihoods
  
  # -- Generate Candidate Points -- #
  if(randomCandidates == T){
    candidates = runif(numCandidates, xmin, xmax)
  } else if(randomCandidates == F){
    candidates = seq(xmin, xmax, length.out = numCandidates)
  } else {
    warning("the randomCandidates argument is not an option; default candidate generation is used.")
  }
  
  # -- Initialize 1st Design Point in D -- #
  if(initialization == 0){
    # initialization == 0: get the point at which f1 and f2 are most different
    xinitind = which.max(abs(f0(candidates) - f1(candidates)))
  } else if(initialization == 1){
    # initialization == 1: get the point at which f1 and f2 are most different
    xinitind = which.max(dnorm(x = candidates, mean = x * beta0, sd = sqrt(var_e)) / 
                                 dnorm(x = candidates, mean = x * beta1, sd = sqrt(var_e)))
  } else{
    warning("the initialization argument is not an option; default initialization is used.")
    xmostdifferentind = which.max(abs(f0(candidates) - f1(candidates)))
  }
  
  D = candidates[xinitind] # x1, first element of set of design points, D
  # also D is used to initialize greedy algorithm
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
    # f_opt = which.min(f_min(candidates, D, k, beta0, beta1, var_e))
    # xnew = candidates[f_opt]
    
    #f_opt = which.min(f_min(candidates, D, k, mean_beta0, mean_beta1, var_e, var_mean))
    f_min_candidates = sapply(candidates, function(x) f_min(x, D, k, beta0, beta1, var_e))
    f_opt = which.min(f_min_candidates)
    xnew = candidates[f_opt]
    
    # Update set of design points (D) and plot new point
    D = c(D,xnew) # add the new point to the set
    candidates = candidates[-f_opt]
    text(xnew, f0(xnew), i, col = 4)
    text(xnew, f1(xnew), i, col = 4)
    points(xnew, 0, col = 2)
  }
  return(D)
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

beta0 = 1 # slope of null model
beta1 = 1 / 2 # slope of alternative model
var_e = 1 # same variance

f0 = function(x) beta0 * x # null regression model
f1 = function(x) beta1 * x # alternative regression model

n = 11
numCandidates = 1000
k = 4 # optimal k was k = 4p
xmin = 0
xmax = 1

X_test = SMED_ms(f0, beta0, f1, beta1, var_e, n, numCandidates, k, xmin, xmax)




#dev.copy(png,'linear_regression.png')
#dev.off()







# testing wasserstein function from "transport" library
x = 0.3
Wasserstein_distance(x, beta0, beta1, var_e)
wasserstein1d(f0(x), f1(x))
# it works!

library(OptimalDesign)
?od.m1
candidates = runif(numCandidates, xmin, xmax)
od.m1(candidates, b, A=NULL, w0=NULL, type="exact", kappa=1e-9,
      tab=NULL, graph=NULL, t.max=120)