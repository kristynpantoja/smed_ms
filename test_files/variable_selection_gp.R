# testing higher dimensions

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


# --- simulations  --- #
numSims = 3

# x_seq, grid over which to generate subsequent functions
xmin = 0; xmax = 1
numx = 1001
x_seq = seq(from = xmin, to = xmax, length.out = numx) # set training points

# input set designs
N = 6; N2 = 15
Ntotal = N + N2

# MMED parameters for testing
l01= c(0.1, 0.1)
# type01 = c(1, 1)
numCandidates = 1001
k = 4
p = 1
nugget = NULL
alpha = 1

# input set 1
set.seed(6)
x_train_ind = sample(1:numx, N)
x_train = x_seq[x_train_ind]


### 1 DIMENSION vs 1 DIMENSION (comparing two different kernels)
type01 = c(1, 4)
# generate matern functions
set.seed(12)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
seed = 1
set.seed(seed)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values


# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)

# run simulations
start_time <- Sys.time()
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  
  # generate mmed
  mmed_gp_list[[i]] = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = N2, 
                                                    k = k, p = p, xmin = xmin, xmax = xmax, 
                                                    nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  
  # mmed inputs' predictions and evaluations
  newpts = mmed_gp_list[[i]]$addD
  newpts_mat[ , i] = newpts
  newpts_ind = mmed_gp_list[[i]]$indices
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  if(i == numSims){
    H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
    H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
    postpredmu0 = H0_pred$pred_mean
    postpredmu1 = H1_pred$pred_mean
    plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
         main = "1D") # plot the function
    points(x = x_train, y = y_train)
    H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
    H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
    lines(x = x_seq, y = H0_predfn$pred_mean)
    lines(x = x_seq, y = H1_predfn$pred_mean)
    err0 = 2 * sqrt(diag(H0_predfn$pred_var))
    err1 = 2 * sqrt(diag(H1_predfn$pred_var))
    add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
    add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
    points(x = newpts, y = truey, col = 2)
    points(x = newpts, y = postpredmu0, col = 3)
    points(x = newpts, y = postpredmu1, col = 4)
    legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
           col = c(1:4), pch = rep(1,4))
  }
}
end_time <- Sys.time()
end_time - start_time

train1sims_1v1 = list("grid" = x_seq,
                     "sim_fns" = y_seq_mat,
                     "x_train" = x_train,
                     "x_train_ind" = x_train_ind,
                     "mmed_gp_list" = mmed_gp_list,
                     "x_mmed" = newpts_mat,
                     "x_mmed_ind" = newpts_ind_mat,
                     "y_mmed" = truey_mat)




##############################
##############################
##############################
##############################
##############################
##############################
##############################
##############################

numx2 = 21
x_seq2 = seq(from = xmin, to = xmax, length.out = numx2)
x_grid = expand.grid(x_seq2, x_seq2)

# using the designs we have will be tricky. just uniformly select points.
set.seed(6)
x_train_ind = sample(1:numx2, N) # figure out a better way to sample! ##############################
x_train = x_grid[x_train_ind + numx2*(1:N), ]
# x_grid[x_train_ind, ] 
# x_grid[x_train_ind + num  x2*(1:N), 1] # same 1st variable values

### 1 DIMENSION vs 2 DIMENSIONS
type01 = c(1, 1)
# generate matern functions (true dim = 2)
set.seed(12)
null_cov = getCov(x_grid, x_grid, type01[2], l01[2])
null_mean = rep(0, numx2^2)
seed = 1
set.seed(seed)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values
sim_index = 3
quilt.plot(x_grid, y_seq_mat[ , sim_index])
plot(x_grid[ , 1], y_seq_mat[ , sim_index])
# generate matern functions (true dim = 1)
set.seed(12)
null_cov = getCov(x_seq2, x_seq2, type01[2], l01[2])
null_mean = rep(0, numx2)
seed = 1
set.seed(seed)
y_seq_mat_1d = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values
# expand to 2 dims
y_seq_mat = matrix(NA, nrow = numx2^2, ncol = numSims)
for(i in 1:numSims){
  expanded = expand.grid(y_seq_mat_1d[ , i], y_seq_mat_1d[ , i])
  y_seq_mat[ , i] = expanded[ , 1]
}
sim_index = 3
quilt.plot(x_grid, y_seq_mat[ , sim_index])
plot(x_grid[ , 1], y_seq_mat[ , sim_index])
order_truef = order(x_grid[,1])
lines(x_grid[order_truef,1], y_seq_mat[order_truef,sim_index])

# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)

# run simulations
start_time <- Sys.time()
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  
  # generate mmed
  mmed_gp_list[[i]] = add_MED_ms_oneatatime_data_gp2(x_train, y_train, type01, l01, subdim = 1, var_e = 1, N2 = N2, 
                                                    k = k, p = p, xmin = xmin, xmax = xmax, 
                                                    nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  
  # mmed inputs' predictions and evaluations
  newpts = mmed_gp_list[[i]]$addD
  newpts_mat[ , i] = newpts
  newpts_ind = mmed_gp_list[[i]]$indices
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  if(i == numSims){
    H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
    H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
    postpredmu0 = H0_pred$pred_mean
    postpredmu1 = H1_pred$pred_mean
    # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
    #      main = "2D") # plot the function
    # points(x = x_train, y = y_train)
    # H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
    # H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
    # lines(x = x_seq, y = H0_predfn$pred_mean)
    # lines(x = x_seq, y = H1_predfn$pred_mean)
    # err0 = 2 * sqrt(diag(H0_predfn$pred_var))
    # err1 = 2 * sqrt(diag(H1_predfn$pred_var))
    # add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
    # add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
    # points(x = newpts, y = truey, col = 2)
    # points(x = newpts, y = postpredmu0, col = 3)
    # points(x = newpts, y = postpredmu1, col = 4)
    # legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
    #        col = c(1:4), pch = rep(1,4))
  }
}
end_time <- Sys.time()
end_time - start_time

train1sims1v2 = list("grid" = x_seq,
                     "sim_fns" = y_seq_mat,
                     "x_train" = x_train,
                     "x_train_ind" = x_train_ind,
                     # save the 3 designs (post input pts) that we're comparing
                     "x_spacefill" = x_spacefill,
                     "x_spacefill_ind" = x_spacefill_ind,
                     "mmed_gp_list" = mmed_gp_list,
                     "x_mmed" = newpts_mat,
                     "x_mmed_ind" = newpts_ind_mat,
                     "y_mmed" = truey_mat)

