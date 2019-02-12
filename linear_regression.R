# y_i = x_i beta + epsilon_i
# epsilon_i ~ N(0, sigma_e^2)
# y_i ~ N(x_i beta, sigma_e^2)

library(transport)

### Helper  Functions ###

# Wasserstein distance d_W betwen two (univariate) normals
Wasserstein_distance = function(x, beta0, beta1, var_e){
  mu1 = x * beta0
  mu2 = x * beta1
  var1 = var_e
  var2 = var_e
  return (mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2)
}

########## incorporate bayes factor, or nah? #######################################################################################
# charge function at design point x
q = function(x, beta0, beta1, var_e) 1.0 / Wasserstein_distance(x, beta0, beta1, var_e)

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min <- Vectorize(function(candidate, D, k, beta0, beta1, var_e) {
  # xnew = an x from candidate set, xall = D, design points
  q(candidate, beta0, beta1, var_e)^k * sum(sapply(D,function(x){(q(x, beta0, beta1, var_e)/abs(x-candidate))^k}))
}, vectorize.args='candidate')



### Iterative Algorithm for SMED for Model Selection ###

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(f0, beta0, f1, beta1, var_e, n = 10, navail = 1000, k = 4, xmin = 0, xmax = 1, include.histogram = F){
  # Function that implements SMED-inspired model selection in 1D
  #  f: parameters of the two models, in the same density family (for now - in future, change Wasserstein)
  #  n: # of pts to select
  #  navail: # of pts to use as candidates
  #  k: power to use for MED. k=4p is default
  #  xmin, xmax: limits on inputs
  #  include.histogram - whether to plot histogram at end
  
  # -- Initialize with mode -- #
  xavail <- runif(navail,xmin,xmax) # or seq(xmin, xmax, length.out = navail)
  # get the point that is highest for both f1 and f2 (has largest ratio of likelihoods, instead, for more medium value? ##########################)
  xmodeind <- which.max(abs(f0(xavail) - f1(xavail)))
  X <- xavail[xmodeind] # this is the "good" design point, x1, first element of set of design points, D
  # also X is used to initialize greedy algorithm
  xavail <- xavail[-xmodeind] # candidate set, for choosing next design point x_{n+1}
  
  # Plot density and (highest lik) points
  curve(f0, from = xmin, to = xmax)
  curve(f1, col = 2, add = TRUE)
  text(X,f0(X),1,col=4)
  text(X,f1(X),1,col=4)
  points(X,0,col=2)
  
#  points(X,0,col=2)
  
  # Sequentially pick rest
  for(i in 2:n) {
    # Find f_opt: minimum of f_min
    f_opt <- which.min(f_min(xavail, X, k, beta0, beta1, var_e))
    xnew <- xavail[f_opt]
    # Update X and plot new point
    X <- c(X,xnew) # add the new point to the set
    xavail <- xavail[-f_opt]
    text(xnew,f0(xnew),i,col=4)
    text(xnew,f1(xnew),i,col=4)
    points(xnew,0,col=2)
  }
  if (include.histogram) hist(X,breaks=20)
  return(X)
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each
# f = c(mu, var)
# testing for evidence of two different slopes

beta0 = 1
beta1 = 1/2
var_e = 1

f0 = function(x) beta0 * x
f1 = function(x) beta1 * x

n = 100
navail = 1000
k = 4
xmin = 0
xmax = 1

X_test = SMED_ms(f0, beta0, f1, beta1, var_e, n, navail, k, xmin, xmax)


# testing wasserstein function from "transport" library
x = 0.3
Wasserstein_distance(x, beta0, beta1, var_e)
wasserstein1d(f0(x), f1(x))






