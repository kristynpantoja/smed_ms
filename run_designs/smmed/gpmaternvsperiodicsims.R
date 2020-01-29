# libraries
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# Computer
# home = "/home/kristyn/Documents/smed_ms"

# Cluster
home = "/scratch/user/kristynp/smed_ms"
output_home = paste(home,"/run_designs/",sep="")

# source files for evaluations

# --- sources to generate MEDs --- #

functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/simulate_seqMED.R", sep = ""))

# --- sources to designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
source(paste(functions_home, "/plot_EPH1.R", sep = ""))
source(paste(functions_home, "/plot_MSE.R", sep = ""))
source(paste(functions_home, "/plot_posterior_variance.R", sep = ""))
source(paste(functions_home, "/postpredyhat_mse_closedform.R", sep = ""))

# --- Sources/Libraries for gaussian process stuff  --- #

source(paste(functions_home, "/med_ms_fns_gp.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))

add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}



## Gaussian vs. Periodic

# x_seq, grid over which to generate subsequent functions
xmin = 0; xmax = 1
numx = 1001
x_seq = seq(from = xmin, to = xmax, length.out = numx) # set training points

# train set designs
N = 6; N2 = 15
Ntotal = 21
# 1. make space-filling design
# space_filling = seq(from = xmin, to = xmax, length.out = Ntotal)
space_filling_ind = c(1, 1 + ((numx - 1)/(Ntotal - 1)) * 1:((numx - 1) / ((numx - 1)/(Ntotal - 1))))
space_filling = x_seq[space_filling_ind]

# train set 1
x_train1_ind = space_filling_ind[1:N]
x_train1 = x_seq[x_train1_ind]
x_spacefill1_ind = space_filling_ind[-c(1:N)]
x_spacefill1 = x_seq[x_spacefill1_ind]
# all.equal(space_filling, c(x_train1, x_spacefill1))

# train set 2
x_train2_ind = space_filling_ind[c(1, 2, 4, 7, 12, 21)]
x_train2 = x_seq[x_train2_ind]
x_spacefill2_ind = space_filling_ind[-c(1, 2, 4, 7, 12, 21)]
x_spacefill2 = x_seq[x_spacefill2_ind]
# all.equal(space_filling, sort(c(x_train2, x_spacefill2)))

# train set 3 (space-filling)
x_train3_ind = c(1, 1 + ((numx - 1)/(N - 1)) * 1:((numx - 1) / ((numx - 1)/(N - 1))))
x_train3 = x_seq[x_train3_ind]
x_spacefill3_ind = space_filling_ind[!(space_filling_ind %in% x_train3_ind)]
x_spacefill3 = x_seq[x_spacefill3_ind]
# all.equal(space_filling, sort(c(x_train3, x_spacefill3)))

# MMED parameters for testing
l01= c(0.1, 0.5); type01 = c(4, 5); numCandidates = 1001; k = 4; p = 1; nugget = NULL; alpha = 1

# generate periodic function
set.seed(13)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
y_seq_periodic = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov)) # the function values

# things that stay the same
numSims = 25
seed = 1
set.seed(seed)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values




## TRAIN SET 1 ##

x_train = x_train1
x_train_ind = x_train1_ind
x_spacefill = x_spacefill1
x_spacefill_ind = x_spacefill1_ind

# these get rewritten for each train set
mmed_gp_list = list()
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
# save log likelihood ratio (posterior predictive distr)
logLR01pred_mmed_vec = rep(NA, numSims)
logLR01pred_sf_vec = rep(NA, numSims)
# save log likelihood ratio (joint distr)
logLR01joint_mmed_vec = rep(NA, numSims)
logLR01joint_sf_vec = rep(NA, numSims)

# see what they look like
# plot(x_seq, y_seq_mat[ , 1], type = "l", ylim = range(y_seq_mat), ylab = "sample function",
#      main = "Matern Kernel") # plot the function
# lines(x_seq, y_seq_mat[ , 2], col = 2); lines(x_seq, y_seq_mat[ , 3], col = 3); lines(x_seq, y_seq_mat[ , 4], col = 4)

# generate mmed and compute RSS0 and RSS1 for mmed and space_fill each
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
  # predictions using each model (using mmed)
  newpts = mmed_gp_list[[i]]$addD
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
  # lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
  # lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
  # err0 = 2 * sqrt(diag(H0_predfn$pred_var))
  # err1 = 2 * sqrt(diag(H1_predfn$pred_var))
  # add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
  # add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
  truey = y_seq[mmed_gp_list[[i]]$indices]
  # points(x = newpts, y = truey, col = 2)
  # points(x = newpts, y = postpredmu0, col = 3)
  # points(x = newpts, y = postpredmu1, col = 4)
  # legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
  #        col = c(1:4), pch = rep(1,4))
  # it looks like the matern is better for modeling this function
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for mmed (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_mmed_vec[i] = likH0 - likH1
  
  # space filling: (note: y_train is different for each function, which is why we need this in the loop)
  # predictions using each model
  newpts = x_spacefill
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  truey = y_seq[x_spacefill_ind]
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for spacefilling (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_sf_vec[i] = likH0 - likH1
  
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
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

RSS01mmed_vec = RSS0mmed_vec/RSS1mmed_vec
RSS01sf_vec = RSS0sf_vec/RSS1sf_vec

train1sims = list("domain" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "mmed_gp_list" = mmed_gp_list,
                  "RSS01mmed_vec" = RSS01mmed_vec,
                  "logLR01pred_mmed_vec" = logLR01pred_mmed_vec,
                  "RSS01sf_vec" = RSS01sf_vec,
                  "logLR01pred_sf_vec" = logLR01pred_sf_vec)
saveRDS(train1sims, paste(output_home,"maternvsperiodic_train1sims.rds", sep = ""))



## TRAIN SET 2 ##

x_train = x_train2
x_train_ind = x_train2_ind
x_spacefill = x_spacefill2
x_spacefill_ind = x_spacefill2_ind

# these get rewritten for each train set
mmed_gp_list = list()
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
# save log likelihood ratio (posterior predictive distr)
logLR01pred_mmed_vec = rep(NA, numSims)
logLR01pred_sf_vec = rep(NA, numSims)
# save log likelihood ratio (joint distr)
logLR01joint_mmed_vec = rep(NA, numSims)
logLR01joint_sf_vec = rep(NA, numSims)

# see what they look like
# plot(x_seq, y_seq_mat[ , 1], type = "l", ylim = range(y_seq_mat), ylab = "sample function",
#      main = "Matern Kernel") # plot the function
# lines(x_seq, y_seq_mat[ , 2], col = 2); lines(x_seq, y_seq_mat[ , 3], col = 3); lines(x_seq, y_seq_mat[ , 4], col = 4)

# generate mmed and compute RSS0 and RSS1 for mmed and space_fill each
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
  # predictions using each model (using mmed)
  newpts = mmed_gp_list[[i]]$addD
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
  # lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
  # lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
  # err0 = 2 * sqrt(diag(H0_predfn$pred_var))
  # err1 = 2 * sqrt(diag(H1_predfn$pred_var))
  # add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
  # add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
  truey = y_seq[mmed_gp_list[[i]]$indices]
  # points(x = newpts, y = truey, col = 2)
  # points(x = newpts, y = postpredmu0, col = 3)
  # points(x = newpts, y = postpredmu1, col = 4)
  # legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
  #        col = c(1:4), pch = rep(1,4))
  # it looks like the matern is better for modeling this function
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for mmed (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_mmed_vec[i] = likH0 - likH1
  
  # space filling: (note: y_train is different for each function, which is why we need this in the loop)
  # predictions using each model
  newpts = x_spacefill
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  truey = y_seq[x_spacefill_ind]
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for spacefilling (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_sf_vec[i] = likH0 - likH1
  
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
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

RSS01mmed_vec = RSS0mmed_vec/RSS1mmed_vec
RSS01sf_vec = RSS0sf_vec/RSS1sf_vec

train2sims = list("domain" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "mmed_gp_list" = mmed_gp_list,
                  "RSS01mmed_vec" = RSS01mmed_vec,
                  "logLR01pred_mmed_vec" = logLR01pred_mmed_vec,
                  "RSS01sf_vec" = RSS01sf_vec,
                  "logLR01pred_sf_vec" = logLR01pred_sf_vec)
saveRDS(train2sims, paste(output_home,"maternvsperiodic_train2sims.rds", sep = ""))



## TRAIN SET 3 ##

x_train = x_train3
x_train_ind = x_train3_ind
x_spacefill = x_spacefill3
x_spacefill_ind = x_spacefill3_ind

# these get rewritten for each train set
mmed_gp_list = list()
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
# save log likelihood ratio (posterior predictive distr)
logLR01pred_mmed_vec = rep(NA, numSims)
logLR01pred_sf_vec = rep(NA, numSims)
# save log likelihood ratio (joint distr)
logLR01joint_mmed_vec = rep(NA, numSims)
logLR01joint_sf_vec = rep(NA, numSims)

# see what they look like
# plot(x_seq, y_seq_mat[ , 1], type = "l", ylim = range(y_seq_mat), ylab = "sample function",
#      main = "Matern Kernel") # plot the function
# lines(x_seq, y_seq_mat[ , 2], col = 2); lines(x_seq, y_seq_mat[ , 3], col = 3); lines(x_seq, y_seq_mat[ , 4], col = 4)

# generate mmed and compute RSS0 and RSS1 for mmed and space_fill each
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
  # predictions using each model (using mmed)
  newpts = mmed_gp_list[[i]]$addD
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
  # lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
  # lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
  # err0 = 2 * sqrt(diag(H0_predfn$pred_var))
  # err1 = 2 * sqrt(diag(H1_predfn$pred_var))
  # add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
  # add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
  truey = y_seq[mmed_gp_list[[i]]$indices]
  # points(x = newpts, y = truey, col = 2)
  # points(x = newpts, y = postpredmu0, col = 3)
  # points(x = newpts, y = postpredmu1, col = 4)
  # legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
  #        col = c(1:4), pch = rep(1,4))
  # it looks like the matern is better for modeling this function
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for mmed (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_mmed_vec[i] = likH0 - likH1
  
  # space filling: (note: y_train is different for each function, which is why we need this in the loop)
  # predictions using each model
  newpts = x_spacefill
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  truey = y_seq[x_spacefill_ind]
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for spacefilling (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_sf_vec[i] = likH0 - likH1
  
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
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

RSS01mmed_vec = RSS0mmed_vec/RSS1mmed_vec
RSS01sf_vec = RSS0sf_vec/RSS1sf_vec

train3sims = list("domain" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "mmed_gp_list" = mmed_gp_list,
                  "RSS01mmed_vec" = RSS01mmed_vec,
                  "logLR01pred_mmed_vec" = logLR01pred_mmed_vec,
                  "RSS01sf_vec" = RSS01sf_vec,
                  "logLR01pred_sf_vec" = logLR01pred_sf_vec)
saveRDS(train3sims, paste(output_home,"maternvsperiodic_train3sims.rds", sep = ""))



## TRAIN SET 4 ##

x_train_mat = matrix(NA, N, numSims)
x_train_ind_mat = matrix(NA, N, numSims)

# these get rewritten for each train set
mmed_gp_list = list()
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
# save log likelihood ratio (posterior predictive distr)
logLR01pred_mmed_vec = rep(NA, numSims)
logLR01pred_sf_vec = rep(NA, numSims)

# generate mmed and compute RSS0 and RSS1 for mmed and space_fill each
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get x_train
  x_train_ind = sample(1:numx, N)
  x_train_ind_mat[ , i] = x_train_ind
  x_train = x_seq[x_train_ind]
  x_train_mat[ , i] = x_train
  x_spacefill_ind = floor(c(1, 1 + ((numx - 1)/(N2 - 1)) * 1:((numx - 1) / ((numx - 1)/(N2 - 1)))))
  x_spacefill = x_seq[x_spacefill_ind]
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  # generate mmed
  mmed_gp_list[[i]] = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = N2, 
                                                    k = k, p = p, xmin = xmin, xmax = xmax, 
                                                    nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  # predictions using each model (using mmed)
  newpts = mmed_gp_list[[i]]$addD
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
  # lines(x = x_seq, y = H0_predfn$pred_mean)
  # lines(x = x_seq, y = H1_predfn$pred_mean)
  # err0 = 2 * sqrt(diag(H0_predfn$pred_var))
  # err1 = 2 * sqrt(diag(H1_predfn$pred_var))
  # add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
  # add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
  truey = y_seq[mmed_gp_list[[i]]$indices]
  # points(x = newpts, y = truey, col = 2)
  # points(x = newpts, y = postpredmu0, col = 3)
  # points(x = newpts, y = postpredmu1, col = 4)
  # legend("bottomright", legend = c("train y", "true y", "H0:gaussian", "H1:periodic"),
  #        col = c(1:4), pch = rep(1,4))
  # it looks like the matern is better for modeling this function
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for mmed (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_mmed_vec[i] = likH0 - likH1
  
  # space filling: (note: y_train is different for each function, which is why we need this in the loop)
  # predictions using each model
  newpts = x_spacefill
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  truey = y_seq[x_spacefill_ind]
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log lik ratio for spacefilling (pred)
  likH0 = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  likH1 = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  logLR01pred_sf_vec[i] = likH0 - likH1
  
  # plot(x_seq, y_seq, type = "l", ylim = range(y_seq, postpredmu0, postpredmu1), ylab = "sample function",
  #    main = "Matern Kernel") # plot the function
  # points(x = x_train, y = y_train)
  # H0_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
  # H1_predfn = getPredDistrSeq(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
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

RSS01mmed_vec = RSS0mmed_vec/RSS1mmed_vec
RSS01sf_vec = RSS0sf_vec/RSS1sf_vec

train4sims = list("domain" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train_mat,
                  "x_train_ind" = x_train_ind_mat,
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "mmed_gp_list" = mmed_gp_list,
                  "RSS01mmed_vec" = RSS01mmed_vec,
                  "logLR01pred_mmed_vec" = logLR01pred_mmed_vec,
                  "RSS01sf_vec" = RSS01sf_vec,
                  "logLR01pred_sf_vec" = logLR01pred_sf_vec)
saveRDS(train4sims, paste(output_home,"maternvsperiodic_train4sims.rds", sep = ""))





