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


# making sure it works: (seems like it is)
xmin = 0; xmax = 1; initN = 5; l = 0.1 # set some parameters
x_train = seq(from = xmin, to = xmax, length.out = initN) # set training points
# get training data: draw from the prior / null distribution with squared exponential covariance fn
null_cov = getCov(x_train, x_train, 1, l)
null_mean = rep(0, initN)
set.seed(125)
y_train = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov))
l= c(0.1, 0.1); type =  c(1, 2); numCandidates = 100; k = 1; p = 2
# test fn
MED_ms_gp_1pt(x_train, y_train, type, l, initN, numCandidates, k, p, xmin, xmax, nugget = NULL)
# check Wasserstein_distance_postpred_gp, setting std dev to be 0
# method 1, with nugget (bc of necessity)
N2 = 11
updateDesign = add_MED_ms_oneatatime_data_gp_old(x_train, y_train, type, l, var_e = 1, N2 = N2, numCandidates, k = k, p = p,
                                             xmin = 0, xmax = 1, nugget = 1e-15, alpha = NULL, buffer = 0, 
                                             genCandidates = 1)
hist(x_train)
hist(updateDesign$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign$addD, y = rep(0, N2), col = 2)
# method 2, leaving out terms in TPE summation corresp to data points
updateDesign2 = add_MED_ms_oneatatime_data_gp(x_train, y_train, type, l, var_e = 1, N2 = N2, numCandidates, k = k, p = p,
                                             xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                             genCandidates = 1)
hist(x_train)
hist(updateDesign2$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign2$addD, y = rep(0, N2), col = 2)
# method 3, like method 2 but with nugget term also
updateDesign3 = add_MED_ms_oneatatime_data_gp(x_train, y_train, type, l, var_e = 1, N2 = N2, numCandidates, k = k, p = p,
                                               xmin = 0, xmax = 1, nugget = 1e-15, alpha = NULL, buffer = 0, 
                                               genCandidates = 1)
hist(x_train)
hist(updateDesign3$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign3$addD, y = rep(0, N2), col = 2)







# looking at q's as x_i -> x in the set
ind = 2 # any of 2:initN are fine
x_ref = x_train[ind]
seq_start = x_train[ind - 1] + (x_train[ind] - x_train[ind - 1]) / 2
x_seq = seq(from = seq_start, to = x_ref, length.out = 500)
# with nugget:
nugget = 1e-15
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]) + diag(rep(nugget, length(x_train))))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]) + diag(rep(nugget, length(x_train))))
w_seq = sapply(x_seq, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l))
q_seq = sapply(x_seq, FUN = function(x) q_data_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l, p=2, alpha = NULL, buffer = 0))
plot(x_seq, w_seq, type = "l")
plot(x_seq, q_seq, type = "l")
which.max(w_seq)
which.min(q_seq)
# without nugget:
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]))
w_seq = sapply(x_seq, FUN = function(x1) Wasserstein_distance_postpred_gp(x1, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l))
q_seq = sapply(x_seq, FUN = function(x) q_data_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l, p=2, alpha = NULL, buffer = 0))
plot(x_seq, w_seq, type = "l")
plot(x_seq, q_seq, type = "l")
which.max(w_seq)
which.min(q_seq)

# nugget does not affect which point is selected next?

# looking at q's as x_i -> x in the set
x_ref = x_train[ind]
seq_start = x_train[ind - 1]
x_seq = seq(from = seq_start, to = x_ref, length.out = 500)
# with nugget:
nugget = 1e-15
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]) + diag(rep(nugget, length(x_train))))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]) + diag(rep(nugget, length(x_train))))
w_seq = sapply(x_seq, FUN = function(x) Wasserstein_distance_postpred_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l))
q_seq = sapply(x_seq, FUN = function(x) q_data_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l, p=2, alpha = NULL, buffer = 0))
plot(x_seq, w_seq, type = "l")
plot(x_seq, q_seq, type = "l")
which.max(w_seq)
which.min(q_seq)
# without nugget:
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]))
w_seq = sapply(x_seq, FUN = function(x1) Wasserstein_distance_postpred_gp(x1, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l))
q_seq = sapply(x_seq, FUN = function(x) q_data_gp(x, Kinv0, Kinv1, x_train, y_train, var_e=1, type, l, p=2, alpha = NULL, buffer = 0))
plot(x_seq, w_seq, type = "l")
plot(x_seq, q_seq, type = "l")
which.max(w_seq)
which.min(q_seq)




# seeing how their RSS's compare for each model
# consider method 2: faster decay -> smoother -> favors squared exponential
hist(x_train)
hist(updateDesign2$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign2$addD, y = rep(0, 11), col = 2)
# from an actual function
# f(x) = e^{-1.4 x} \text{cos}(7 \pi x / 2), 0 \leq x \leq 1
f1 = function(x) exp(-5* x) * cos(7 * pi * x / 2)
xmin = 0; xmax = 1; N = 5; l = 0.1 # set some parameters
x_train = seq(from = xmin, to = xmax, length.out = N) # set training points
# get training data: from function f1
null_cov = getCov(x_train, x_train, 1, l)
null_mean = rep(0, N)
y_train = sapply(x_train, FUN = f1)
l= c(0.1, 0.1); type =  c(1, 4); numCandidates = 100; k = 4; p = 2
curve(f1)
points(x = x_train, y = y_train)
# predictions using each model
newpts = updateDesign2$addD
k0 = t(as.matrix(getCov(newpts, x_train, type[1], l[1])))
k1 = t(as.matrix(getCov(newpts, x_train, type[2], l[2])))
# (without nugget)
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]))
postpredmu0 = t(k0) %*% Kinv0 %*% y_train
postpredmu1 = t(k1) %*% Kinv1 %*% y_train
curve(f1)
points(x = x_train, y = y_train)
points(x = newpts, y = postpredmu0, col = 3)
points(x = newpts, y = postpredmu1, col = 4)
# it looks like the matern is better for modeling this function
truey = sapply(newpts, FUN = f1)
RSS0 = 1/length(newpts) * sum((postpredmu0 - truey)^2)  
RSS1 = 1/length(newpts) * sum((postpredmu1 - truey)^2)  

# from another function: slower decay -> more erratic/wiggly -> favors matern
# f(x) = e^{-1.4 x} \text{cos}(7 \pi x / 2), 0 \leq x \leq 1
f1 = function(x) exp(-1.2* x) * cos(7 * pi * x / 2)
xmin = 0; xmax = 1; N = 5; l = 0.1 # set some parameters
x_train = seq(from = xmin, to = xmax, length.out = N) # set training points
# get training data: from function f1
null_cov = getCov(x_train, x_train, 1, l)
null_mean = rep(0, N)
y_train = sapply(x_train, FUN = f1)
l= c(0.1, 0.1); type =  c(1, 4); numCandidates = 100; k = 4; p = 2
curve(f1)
points(x = x_train, y = y_train)
# predictions using each model
newpts = updateDesign$addD
k0 = t(as.matrix(getCov(newpts, x_train, type[1], l[1])))
k1 = t(as.matrix(getCov(newpts, x_train, type[2], l[2])))
# (without nugget)
Kinv0 = solve(getCov(x_train, x_train, type[1], l[1]))
Kinv1 = solve(getCov(x_train, x_train, type[2], l[2]))
postpredmu0 = t(k0) %*% Kinv0 %*% y_train
postpredmu1 = t(k1) %*% Kinv1 %*% y_train
curve(f1)
points(x = x_train, y = y_train)
points(x = newpts, y = postpredmu0, col = 3)
points(x = newpts, y = postpredmu1, col = 4)
# it looks like the matern is better for modeling this function
truey = sapply(newpts, FUN = f1)
RSS0 = 1/length(newpts) * sum((postpredmu0 - truey)^2)  
RSS1 = 1/length(newpts) * sum((postpredmu1 - truey)^2)  

# compare to space-filling?
# length(newpts)
# spacefill = seq(from = xmin, to = xmax, length.out = length(newpts) + 3)
# x_train %in% spacefill
# spacefill[!(spacefill %in% x_train)]











# get a function from the null distribution instead
# from periodic kernel, get a function
set.seed(123)
xmin = 0; xmax = 1; l = 1 # set some parameters
numx = 1001
x_seq = seq(from = xmin, to = xmax, length.out = numx) # set training points
# get training data: draw from the prior / null distribution with squared exponential covariance fn
null_cov = getCov(x_seq, x_seq, 5, l)
null_mean = rep(0, numx)
set.seed(124)
y_seq = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov)) # the function values
plot(x_seq, y_seq, type = "l") # plot the function
# get training set
initN = 5 # must be larger than 3
interval_length = ceiling((numx - 1) / (initN - 1))
spacefilling_ind = c(1, 1 + interval_length * 1:(initN - 2), numx)
x_train = x_seq[spacefilling_ind]
y_train = y_seq[spacefilling_ind]
points(x_train, y_train) # plot the training set

# mmed for gp:
N2 = initN - 1 # 4 * initN - 3 also might work?
interval_length = ceiling((numx - 1) / (initN + N2 - 1))
update_spacefilling_ind = c(1, 1 + interval_length * 1:(initN + N2 - 2), numx)
diff(update_spacefilling_ind)
update_spacefilling = update_spacefilling_ind[!(update_spacefilling_ind %in% spacefilling_ind)]
points(x_seq[update_spacefilling], y = rep(-1, length(update_spacefilling)), col = 3)

type = c(1,5); l = c(1, 1); k = 1
updateDesign = add_MED_ms_oneatatime_data_gp(x_train, y_train, type, l, var_e = 1, N2 = N2, numx, k = 4, p = 2,
                                             xmin = 0, xmax = 1, nugget = 1e-15, alpha = NULL, buffer = 0, 
                                             genCandidates = 1)
hist(updateDesign$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign$addD, y = rep(-1.5, N2), col = 2)
# method 2, leaving out terms in TPE summation corresp to data points
updateDesign2 = add_MED_ms_oneatatime_data_gp2(x_train, y_train, type, l, var_e = 1, N2 = N2, numx, k = 4, p = 2,
                                               xmin = 0, xmax = 1, nugget = NULL, alpha = NULL, buffer = 0, 
                                               genCandidates = 1)
hist(x_train)
hist(updateDesign2$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign2$addD, y = rep(0, N2), col = 2)
# method 3, like method 2 but with nugget term also
updateDesign3 = add_MED_ms_oneatatime_data_gp2(x_train, y_train, type, l, var_e = 1, N2 = N2, numx, k = 4, p = 2,
                                               xmin = 0, xmax = 1, nugget = 1e-15, alpha = NULL, buffer = 0, 
                                               genCandidates = 1)
hist(x_train)
hist(updateDesign3$updatedD)
plot(x = x_train, y = y_train)
points(x = updateDesign3$addD, y = rep(0, N2), col = 2)





# it's probbaly easier to make the big space filling design first,
# then take out the elements for the train set
# like literally





