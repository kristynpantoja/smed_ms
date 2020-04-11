# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/generate_MMEDgp_oneatatime.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(fields)

add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}

# some settings
# for type:
#   1 = gaussian kernel
#   2 = matern kernel (from a package, where length-scale parameter is assumed 1, so I don't use it)
#   3 = exponential kernel
#   4 = matern kernel (with nu = 3 / 2, derived)
#   5 = periodic kernel (with fixed period)
# for all of the above, the only parameter I let to be specified is the length-scale parameter, l
true_type = 4 # matern, i.e. misspecified kernel ########################################################################## CAN ADJUST
true_l = 0.1 ############################################################################################################## CAN ADJUST
hypothesized_type = 1 # gaussian, i.e. misspecified kernel ################################################################ CAN ADJUST
# hypothesized_type = 4 # matern, i.e. correctly specified kernel
hypothesized_l = 0.1 # set length-scale parameter ######################################################################### CAN ADJUST
# design space
xmin = 0
xmax = 1
# number of training points (observed data) to use
N = 6
seed = 12 ################################################################################################################# CAN ADJUST

# if you want to see the parameters used to generate true functions in my experiments:
# hypothesis_type
#   1 = in gaussian vs. matern, the true function is generated from matern kernel
#       with length-scale parameter l = 0.1
#   2 = in matern vs. periodic, the true function is generated from periodic kernel
#       with l = 0.5
hypothesis_type = 0 ####################################################################################################### CAN ADJUST
if(hypothesis_type == 1){
  hypothesized_type = 1 # gaussian
  hypothesized_l = 0.1 # set length-scale parameter
  true_type = 4 # matern
  true_l = 0.1
  seed = 12
}
if(hypothesis_type == 2){
  hypothesized_type = 4 # gaussian
  hypothesized_l = 0.1 # set length-scale parameter
  true_type = 5 # periodic
  true_l = 0.5
  seed = 13
}

# generate function
numx = 1001 # number of candidate points (don't change)
x_seq = seq(from = xmin, to = xmax, length.out = numx)

# generate a function from specified true_type with specified hyperparameter l
null_cov = getCov(x_seq, x_seq, true_type, true_l)
null_mean = rep(0, numx)
# this gives the true function
set.seed(seed)
y_seq = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov)) # the function values
# get training data (x)
# train_set
# 1 = design points in one region of design space
# 2 = design points increasingly spread out
# 3 = space_filling
# 4 = (uniformly) randomly-selected points

# make space-filling design, to be used for selecting input points and comparing to train sets
# space_filling = seq(from = xmin, to = xmax, length.out = Ntotal)
N = 6; N2 = 15
Ntotal = 21
space_filling_ind = c(1, 1 + ((numx - 1)/(Ntotal - 1)) * 1:((numx - 1) / ((numx - 1)/(Ntotal - 1))))
space_filling = x_seq[space_filling_ind]

train_set = 4 ######################################################################################################## CAN ADJUST
train_seed = 123
if(train_set == 1){
  x_train_ind = space_filling_ind[1:N]
  x_train = x_seq[x_train_ind]
}
if(train_set == 2){ # only works when N = 6 in this implementation
  Ntotal = 21
  space_filling_ind = c(1, 1 + ((numx - 1)/(Ntotal - 1)) * 1:((numx - 1) / ((numx - 1)/(Ntotal - 1))))
  x_train_ind = space_filling_ind[c(1, 2, 4, 7, 12, 21)]
  x_train = x_seq[x_train_ind]
}
if(train_set == 3){
  x_train_ind = c(1, 1 + ((numx - 1)/(N - 1)) * 1:((numx - 1) / ((numx - 1)/(N - 1))))
  x_train = x_seq[x_train_ind]
}
if(train_set == 4){
  set.seed(train_seed)
  x_train_ind = sample(1:numx, N)
  x_train = x_seq[x_train_ind]
}
# get training data (y)
y_train = y_seq[x_train_ind] # get corresponding y-values, for training

# get predicted curve, i.e. posterior mean and also its covariance, from hypothesized_type,
#   which may be misspecified.
y_pred = getGPPredictive(x_seq, x_train, y_train, hypothesized_type, hypothesized_l, nugget = NULL)

# plot true function, y_seq, and the predicted function, y_pred
err = 2 * sqrt(diag(y_pred$pred_var)) # for pointwise error
ylim1 = range(y_seq, min(y_pred$pred_mean - err), max(y_pred$pred_mean + err))
plot(x_seq, y_seq, type = "l", ylim = ylim1, ylab = "f(x)", main = "True Function and Reconstructed Function") # plot the function
points(x_train, y_train)
add_errorbands(x_seq, y_pred$pred_mean, err, rgb(1, 0, 0, 0.1))
lines(x_seq, y_pred$pred_mean, col = 2, lty = 3)
legend("bottomright", legend = c("true f(x)", "reconstructed f(x)"), col = c(1, 2), lty = c(1, 3), lwd = c(1, 1))












