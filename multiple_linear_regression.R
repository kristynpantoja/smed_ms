### --- Linear Regression, More General --- ###

# Null Model: y_i = x_i' beta0 + epsilon0_i
# vs. 
# Alternative Model: y_i = x_i' beta1 + epsilon1_i

# where, in both models,
#  epsilon_i ~ N(0, sigma_e^2), sigma_e = {sigma_e0, sigma_e1}
#  y_i ~ N(x_i' beta, sigma_e^2), beta = {beta0, beta1}




### Helper  Functions ###

# Wasserstein distance betwen two (univariate) normals at point x
Wasserstein_distance = function(x, beta0, beta1, var0, var1){
  mu1 = as.numeric(crossprod(x, beta0))
  mu2 = as.numeric(crossprod(x, beta1))
  var1 = var0
  var2 = var1
  return (mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2)
}

# charge function at design point x
q = function(x, beta0, beta1, var0, var1) 1.0 / Wasserstein_distance(x, beta0, beta1, var0, var1)

# minimizing criterion for greedy algorithm - calculate this for each x_i in candidate set,
# select the one with the smallest value of f_min to add to current design set D
f_min <- Vectorize(function(candidate, D, k, beta0, beta1, var0, var1) { # DOESN"T WORK. WILL FIX SOON.
  # xnew = an x from candidate set, xall = D, design points
  q(candidate, beta0, beta1, var0, var1)^k * 
    sum(sapply(D, function(x) (q(x, beta0, beta1, var0, var1) / abs(x - candidate))^k))
}, vectorize.args='candidate')

### Iterative Algorithm for SMED for Model Selection ###
# Function that implements SMED-inspired model selection of Joseph, Dasgupta, Tuo, Wu

# f is normal, so need to specify f0 = {beta0}, f1 = {beta1}, var_e
# since these determine both f0 and f1's normal parameters 
SMED_ms = function(f0, beta0, f1, beta1, var0, var1, n = 10, numCandidates = 1000, k = 4, xmin = NULL, xmax = NULL, 
                   randomCandidates = T, initialization = 0){
  #  f0, beta0, var0 : regression line, slope, and error variance of null model
  #  f1, beta1, var1 : regression line, slope, and error variance of alternative model
  #  n : number of design points to select (for set of design points, D)
  #  numCandidates : # of points to use as candidates (set from which design points are selected)
  #  k : power to use for MED. k = 4p is default
  #  xmin, xmax : limits on input vectors
  # randomCandidates : generate candidate points using uniform distribution by default.
  #                    if FALSE, choose equally-spaced out sequence of candidates
  # initialization : method of initializing first design point of D
  #                  0 = candidate that outputs the largest absolute difference between f0 and f1
  #                  1 = candidate that has the largest ratio of likelihoods
  
  p = length(beta0)
  
  # -- Check xmin and xmax -- #
  if(is.null(xmin) & is.null(xmax)){
    xmin = rep(0, p)
    xmax = rep(1, p)
  } else{
    if(length(xmin) != p) stop("xmin, the lower bounds of the input, does not match dimension of input")
    if(length(xmax) != p) stop("xmax, the upper bounds of the input, does not match dimension of input")
  }
  
  # -- Generate Candidate Points -- #
  candidates = matrix(rep(NA, numCandidates * p), nrow = numCandidates, ncol = p)
  if(randomCandidates == T){
    for(i in 1:p){
      candidates[ , i] = runif(numCandidates, xmin[i], xmax[i])
    }
  } else if(randomCandidates == F){
    for(i in 1:p){
      candidates[ , i] = seq(xmin[i], xmax[i], length.out = numCandidates)
    }
  } else {
    warning("the randomCandidates argument is not an option; default candidate generation is used.")
  }
  
  # -- Initialize 1st Design Point in D -- #
  if(initialization == 0){
    # initialization == 0: get the point at which f1 and f2 are most different
    xinitind = which.max(abs(apply(candidates, 1, f0) - apply(candidates, 1, f1)))
  } else if(initialization == 1){
    # initialization == 1: get the point at which f1 and f2 are most different
    xinitind = which.max(apply(candidates, 1, 
                               function(x) dnorm(x, mean = as.numeric(crossprod(x, beta0)), sd = sqrt(var0)) /
                                 dnorm(x, mean = as.numeric(crossprod(x, beta1)), sd = sqrt(var1))))
  } else{
    warning("the initialization argument is not an option; default initialization is used.")
    xinitind = which.max(abs(apply(candidates, 1, f0) - apply(candidates, 1, f1)))
  }
  D = candidates[xinitind, ] # x1, first element of set of design points, D
  # also D is used to initialize greedy algorithm
  candidates <- candidates[-xinitind, ] # candidate set, for choosing next design point x_{n+1}
  
  # Plot density and (highest lik) points
  if(p == 1){
    curve(f0, from = xmin[1], to = xmax[1])
    curve(f1, col = 2, add = TRUE)
    text(D, f0(D), 1, col=4)
    text(D, f1(D), 1, col=4)
    points(D, 0, col=2)
  } else if(p == 2){
    ### SOMETHING HERE ###############################################################
  }
  
  # Sequentially pick rest
  for(i in 2:n){
    # Find f_opt: minimum of f_min
    f_opt = which.min(f_min(candidates, D, k, beta0, beta1, var0, var1)) # some issue with the vectorize part of f_min
    xnew = candidates[f_opt, ] # subscript out of bounds for some reason?
    # Update set of design points (D) and plot new point
    D = cbind(D,xnew) # add the new point to the set
    candidates = candidates[-f_opt, ]
    if(p == 1){
      text(xnew, f0(xnew), i, col = 4)
      text(xnew, f1(xnew), i, col = 4)
      points(xnew, 0, col = 2)
    } else if(p == 2){
      ### SOMETHING HERE ###############################################################
    }
  }
  return(D)
}

# criterion (1 / d_W(f0, f1; x_i) * 1 / d_W(f0, f1; x_j)) / d(x_i, x_j)
# here, we assume f's are normal, so we specify the mean and variance of each

p = 2

beta0 = rep(1, p) # vectors of parameters of null model
beta1 = rep(1 / 2, p) # vectors of parameters of alternative model
var0 = rep(1, 1)
var1 = var0 # same variance

f0 = function(x) as.numeric(crossprod(beta0, x)) # null regression model
f1 = function(x) as.numeric(crossprod(beta1, x)) # alternative regression model

n = 10
numCandidates = 500
k = 4 * p
xmin = rep(0, p)
xmax = rep(1, p)

X_test = SMED_ms(f0, beta0, f1, beta1, var0, var1, n, numCandidates, k, xmin, xmax)


# testing wasserstein function from "transport" library
x = 0.3
Wasserstein_distance(x, beta0, beta1, var_e)
wasserstein1d(f0(x), f1(x))
# it works!