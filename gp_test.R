# --- Working Directory --- #
home = "/Users/kristyn/Documents/research/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/med_ms_fns_gp.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(reshape2)
library(ggplot2)
library(mvtnorm)
library(RandomFieldsUtils)

# library(transport)
# library(mined)


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################   Gaussian Process Model Selection RETRY   ####################
####################################################################################
####################################################################################
####################################################################################
####################################################################################




# Suppose the function is from a GP with squared exponential, want to use predictive posterior
# to get more design points
# Start off with training data (generated from a draw from null distribution with squared exponential k)
xmin = 0; xmax = 1; N = 11; l = 0.1 # set some parameters
x_train = seq(from = xmin, to = xmax, length.out = N) # set training points
# get training data: draw from the prior / null distribution with squared exponential covariance fn
null_cov = getCov(x_train, x_train, 1, l)
null_mean = rep(0, N)
set.seed(124)
y_train = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov))
plot(x_train, y_train) # this is our data, which we pretend we started off with.

# what does predictive distribution look like, for a new point x_star?
x_star = 0.25
# getPredDistr = function(x_star, x_train, y_train, l, C_fn_type){
#   k_star = t(as.matrix(C_fn_1d(x_star, x_train, l, C_fn_type)))
#   K_obs = C_fn_1d(x_train, x_train, l, C_fn_type)
#   pred_mean = t(k_star) %*% solve(K_obs, y_train)
#   pred_cov = C_fn_1d(x_star, x_star, l, C_fn_type) - t(k_star) %*% solve(K_obs, k_star)
#   return(list("pred_mean" = pred_mean, "pred_cov" = pred_cov))
# }
pred_dist_1 = getPredDistr(x_star, x_train, y_train, 1, l, nugget = NULL)
pred_dist_2 = getPredDistr(x_star, x_train, y_train, 2, l, nugget = NULL)
pred_dist_3 = getPredDistr(x_star, x_train, y_train, 3, l, nugget = NULL)
pred_dist_4 = getPredDistr(x_star, x_train, y_train, 4, l, nugget = NULL)
plot(x_train, y_train) # this is our data, which we pretend we started off with.
points(x_star, pred_dist_1$pred_mean, col = 2)
points(x_star, pred_dist_2$pred_mean, col = 2)
points(x_star, pred_dist_3$pred_mean, col = 2)
points(x_star, pred_dist_4$pred_mean, col = 2)


# suppose H0: C_fn_type = 1 = squared exponential
#         H1: C_fn_type = 3 = matern with v = 3/2

# Step 1: get candidate set
numCandidates = 100
x_cand = seq(from = xmin, to = xmax, length.out = numCandidates)
#plot(x_cand, pch = ".")

# Step 2: pick a new design point from the candidate set
# 



###########
# test fn #
###########

home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/med_ms_fns_gp.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/med_ms_fns_gp.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))

xmin = 0; xmax = 1; N = 5; l = 0.1 # set some parameters
x_train = seq(from = xmin, to = xmax, length.out = N) # set training points
# get training data: draw from the prior / null distribution with squared exponential covariance fn
null_cov = getCov(x_train, x_train, 1, l)
null_mean = rep(0, N)
set.seed(125)
y_train = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov))
l= c(0.1, 0.1); type =  c(1, 4); numCandidates = 100; k = 4; p = 2
# test fn
MED_ms_gp_1pt(x_train, y_train, type, l, N, numCandidates, k, p, xmin, xmax, nugget = NULL)
# check Wasserstein_distance_postpred_gp, setting std dev to be 0
updateDesign = add_MED_ms_oneatatime_data_gp(x_train, y_train, type, l, var_e = 1, N2 = 11, numCandidates, k = 4, p = 2,
                                             xmin = 0, xmax = 1, nugget = 1e-15, alpha = NULL, buffer = 0, 
                                             genCandidates = 1)
hist(x_train)
hist(updateDesign$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign$addD, y = rep(0, 11), col = 2)

updateDesign2 = add_MED_ms_oneatatime_data_gp2(x_train, y_train, type, l, var_e = 1, N2 = 11, numCandidates, k = 4, p = 2,
                                             xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                             genCandidates = 1)
hist(x_train)
hist(updateDesign2$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign2$addD, y = rep(0, 11), col = 2)


# leave out second summation, with x_i in dataset
# consider doing the same for SMMED (for polynomial models) to see if the results are the same??
#   I don't see why they would be, since then the distances between candidates x and points x_i in initD
#   are not being considered
# calculate RSS for y-hat from K0 and y-hat from K1

