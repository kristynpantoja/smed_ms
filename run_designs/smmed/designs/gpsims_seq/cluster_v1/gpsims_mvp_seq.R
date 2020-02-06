##################################################
## --- Gaussian vs. Periodic (Sequentially) --- ##
##################################################

# --- libraries --- #

# libraries
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# --- workspaces --- #

# Computer
# home = "/home/kristyn/Documents/smed_ms"
# output_home = paste(home, "/", sep = "")

# Cluster
home = "/scratch/user/kristynp/smed_ms"
output_home = paste(home,"/run_designs_v1/",sep="")

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

# --- sources/functions for gaussian process stuff  --- #

source(paste(functions_home, "/med_ms_fns_gp.R", sep = ""))
source(paste(functions_home, "/update_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))

add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}

# --- simulations  --- #
numSims = 25

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

# generate periodic functions
set.seed(13)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
seed = 1
set.seed(seed)
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values

# sequential stuff
numSteps = 3
stepN = 5 # get 5 new points at each step



## INPUT SET 1 ##

# input set 1
x_train = x_train1
x_train_ind = x_train1_ind

# space-filling for input set 1
x_spacefill = x_spacefill1
x_spacefill_ind = x_spacefill1_ind

# these get rewritten for each train set
# uniform
x_uniform_mat = matrix(NA, N2, numSims)
x_uniform_ind_mat = matrix(NA, N2, numSims)
# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# RSS
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
RSS0unif_vec = rep(NA, numSims)
RSS1unif_vec = rep(NA, numSims)
# log predictive density
logPredEvid0mmed_vec = rep(NA, numSims)
logPredEvid1mmed_vec = rep(NA, numSims)
logPredEvid0sf_vec = rep(NA, numSims)
logPredEvid1sf_vec = rep(NA, numSims)
logPredEvid0unif_vec = rep(NA, numSims)
logPredEvid1unif_vec = rep(NA, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)
logJointEvid0unif_vec = rep(NA, numSims)
logJointEvid1unif_vec = rep(NA, numSims)

# run simulations
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  # get x_uniform
  x_uniform_ind = sample(1:numx, N2)
  x_uniform_ind_mat[ , i] = x_uniform_ind
  x_uniform = x_seq[x_uniform_ind]
  x_uniform_mat[ , i] = x_uniform
  
  # generate mmed
  # step 1:
  step1_mmed = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step1_newpts = step1_mmed$addD
  
  x_train_new = c(x_train, step1_mmed$addD)
  x_train_ind_new = c(x_train_ind, step1_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 2:
  step2_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step2_newpts = step2_mmed$addD
  
  x_train_new = c(x_train_new, step2_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step2_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 3:
  step3_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step3_newpts = step3_mmed$addD
  
  x_train_new = c(x_train_new, step3_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step3_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  
  mmed_gp_list[[i]] = list(step1_mmed, step2_mmed, step3_mmed)
  
  # mmed inputs' predictions and evaluations
  newpts = c(step1_newpts, step2_newpts, step3_newpts)
  newpts_mat[ , i] = newpts
  newpts_ind = c(step1_mmed$indices, step2_mmed$indices, step3_mmed$indices)
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for mmed
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for mmed
  logJointEvid0mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # space-filling inputs' predictions and evaluations
  # (note: y_train is different for each function, which is why we need this in the loop)
  newpts = x_spacefill
  truey = y_seq[x_spacefill_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for spacefilling
  logPredEvid0sf_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1sf_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for spacefilling
  logJointEvid0sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # uniform inputs' predictions and evaluations
  newpts = x_uniform
  truey = y_seq[x_uniform_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for uniform
  RSS0unif_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1unif_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for uniform
  logPredEvid0unif_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1unif_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for uniform
  logJointEvid0unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
}

train1sims = list("grid" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  # save the 3 designs (post input pts) that we're comparing
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "x_uniform" = x_uniform_mat,
                  "x_uniform_ind" = x_uniform_ind_mat,
                  "mmed_gp_list" = mmed_gp_list,
                  "x_mmed" = newpts_mat,
                  "x_mmed_ind" = newpts_ind_mat,
                  "y_mmed" = truey_mat,
                  # save evaluations for each of them also
                  "RSS0mmed_vec" = RSS0mmed_vec,
                  "RSS1mmed_vec" = RSS1mmed_vec,
                  "RSS0sf_vec" = RSS0sf_vec,
                  "RSS1sf_vec" = RSS1sf_vec,
                  "RSS0unif_vec" = RSS0unif_vec,
                  "RSS1unif_vec" = RSS1unif_vec,
                  "RSS01mmed" = RSS0mmed_vec / RSS1mmed_vec,
                  "RSS01sf" = RSS0sf_vec / RSS1sf_vec,
                  "RSS01unif" = RSS0unif_vec / RSS1unif_vec,
                  #
                  "logPredEvid0mmed_vec" = logPredEvid0mmed_vec,
                  "logPredEvid1mmed_vec" = logPredEvid1mmed_vec,
                  "logPredEvid0sf_vec" = logPredEvid0sf_vec,
                  "logPredEvid1sf_vec" = logPredEvid1sf_vec,
                  "logPredEvid0unif_vec" = logPredEvid0unif_vec,
                  "logPredEvid1unif_vec" = logPredEvid1unif_vec,
                  "logPredEvid01mmed" = logPredEvid0mmed_vec - logPredEvid1mmed_vec,
                  "logPredEvid01sf" = logPredEvid0sf_vec - logPredEvid1sf_vec,
                  "logPredEvid01unif" = logPredEvid0unif_vec - logPredEvid1unif_vec,
                  #
                  "logJointEvid0mmed_vec" = logJointEvid0mmed_vec,
                  "logJointEvid1mmed_vec" = logJointEvid1mmed_vec,
                  "logJointEvid0sf_vec" = logJointEvid0sf_vec,
                  "logJointEvid1sf_vec" = logJointEvid1sf_vec,
                  "logJointEvid0unif_vec" = logJointEvid0unif_vec,
                  "logJointEvid1unif_vec" = logJointEvid1unif_vec,
                  "logJointEvid01mmed" = logJointEvid0mmed_vec - logJointEvid1mmed_vec,
                  "logJointEvid01sf" = logJointEvid0sf_vec - logJointEvid1sf_vec,
                  "logJointEvid01unif" = logJointEvid0unif_vec - logJointEvid1unif_vec)
saveRDS(train1sims, paste(output_home,"mvp_seq_train1sims.rds", sep = ""))


## INPUT SET 2 ##

# input set 2
x_train = x_train2
x_train_ind = x_train2_ind

# space-filling for input set 2
x_spacefill = x_spacefill2
x_spacefill_ind = x_spacefill2_ind

# these get rewritten for each train set
# uniform
x_uniform_mat = matrix(NA, N2, numSims)
x_uniform_ind_mat = matrix(NA, N2, numSims)
# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# RSS
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
RSS0unif_vec = rep(NA, numSims)
RSS1unif_vec = rep(NA, numSims)
# log predictive density
logPredEvid0mmed_vec = rep(NA, numSims)
logPredEvid1mmed_vec = rep(NA, numSims)
logPredEvid0sf_vec = rep(NA, numSims)
logPredEvid1sf_vec = rep(NA, numSims)
logPredEvid0unif_vec = rep(NA, numSims)
logPredEvid1unif_vec = rep(NA, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)
logJointEvid0unif_vec = rep(NA, numSims)
logJointEvid1unif_vec = rep(NA, numSims)

# run simulations
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  # get x_uniform
  x_uniform_ind = sample(1:numx, N2)
  x_uniform_ind_mat[ , i] = x_uniform_ind
  x_uniform = x_seq[x_uniform_ind]
  x_uniform_mat[ , i] = x_uniform
  
  # generate mmed
  # step 1:
  step1_mmed = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step1_newpts = step1_mmed$addD
  
  x_train_new = c(x_train, step1_mmed$addD)
  x_train_ind_new = c(x_train_ind, step1_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 2:
  step2_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step2_newpts = step2_mmed$addD
  
  x_train_new = c(x_train_new, step2_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step2_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 3:
  step3_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step3_newpts = step3_mmed$addD
  
  x_train_new = c(x_train_new, step3_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step3_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  
  mmed_gp_list[[i]] = list(step1_mmed, step2_mmed, step3_mmed)
  
  # mmed inputs' predictions and evaluations
  newpts = c(step1_newpts, step2_newpts, step3_newpts)
  newpts_mat[ , i] = newpts
  newpts_ind = c(step1_mmed$indices, step2_mmed$indices, step3_mmed$indices)
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for mmed
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for mmed
  logJointEvid0mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # space-filling inputs' predictions and evaluations
  # (note: y_train is different for each function, which is why we need this in the loop)
  newpts = x_spacefill
  truey = y_seq[x_spacefill_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for spacefilling
  logPredEvid0sf_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1sf_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for spacefilling
  logJointEvid0sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # uniform inputs' predictions and evaluations
  newpts = x_uniform
  truey = y_seq[x_uniform_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for uniform
  RSS0unif_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1unif_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for uniform
  logPredEvid0unif_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1unif_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for uniform
  logJointEvid0unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
}

train2sims = list("grid" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  # save the 3 designs (post input pts) that we're comparing
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "x_uniform" = x_uniform_mat,
                  "x_uniform_ind" = x_uniform_ind_mat,
                  "mmed_gp_list" = mmed_gp_list,
                  "x_mmed" = newpts_mat,
                  "x_mmed_ind" = newpts_ind_mat,
                  "y_mmed" = truey_mat,
                  # save evaluations for each of them also
                  "RSS0mmed_vec" = RSS0mmed_vec,
                  "RSS1mmed_vec" = RSS1mmed_vec,
                  "RSS0sf_vec" = RSS0sf_vec,
                  "RSS1sf_vec" = RSS1sf_vec,
                  "RSS0unif_vec" = RSS0unif_vec,
                  "RSS1unif_vec" = RSS1unif_vec,
                  "RSS01mmed" = RSS0mmed_vec / RSS1mmed_vec,
                  "RSS01sf" = RSS0sf_vec / RSS1sf_vec,
                  "RSS01unif" = RSS0unif_vec / RSS1unif_vec,
                  #
                  "logPredEvid0mmed_vec" = logPredEvid0mmed_vec,
                  "logPredEvid1mmed_vec" = logPredEvid1mmed_vec,
                  "logPredEvid0sf_vec" = logPredEvid0sf_vec,
                  "logPredEvid1sf_vec" = logPredEvid1sf_vec,
                  "logPredEvid0unif_vec" = logPredEvid0unif_vec,
                  "logPredEvid1unif_vec" = logPredEvid1unif_vec,
                  "logPredEvid01mmed" = logPredEvid0mmed_vec - logPredEvid1mmed_vec,
                  "logPredEvid01sf" = logPredEvid0sf_vec - logPredEvid1sf_vec,
                  "logPredEvid01unif" = logPredEvid0unif_vec - logPredEvid1unif_vec,
                  #
                  "logJointEvid0mmed_vec" = logJointEvid0mmed_vec,
                  "logJointEvid1mmed_vec" = logJointEvid1mmed_vec,
                  "logJointEvid0sf_vec" = logJointEvid0sf_vec,
                  "logJointEvid1sf_vec" = logJointEvid1sf_vec,
                  "logJointEvid0unif_vec" = logJointEvid0unif_vec,
                  "logJointEvid1unif_vec" = logJointEvid1unif_vec,
                  "logJointEvid01mmed" = logJointEvid0mmed_vec - logJointEvid1mmed_vec,
                  "logJointEvid01sf" = logJointEvid0sf_vec - logJointEvid1sf_vec,
                  "logJointEvid01unif" = logJointEvid0unif_vec - logJointEvid1unif_vec)
saveRDS(train2sims, paste(output_home,"mvp_seq_train2sims.rds", sep = ""))


## INPUT SET 3 ##

# input set 3
x_train = x_train3
x_train_ind = x_train3_ind

# space-filling for input set 3
x_spacefill = x_spacefill3
x_spacefill_ind = x_spacefill3_ind

# these get rewritten for each train set
# uniform
x_uniform_mat = matrix(NA, N2, numSims)
x_uniform_ind_mat = matrix(NA, N2, numSims)
# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# RSS
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
RSS0unif_vec = rep(NA, numSims)
RSS1unif_vec = rep(NA, numSims)
# log predictive density
logPredEvid0mmed_vec = rep(NA, numSims)
logPredEvid1mmed_vec = rep(NA, numSims)
logPredEvid0sf_vec = rep(NA, numSims)
logPredEvid1sf_vec = rep(NA, numSims)
logPredEvid0unif_vec = rep(NA, numSims)
logPredEvid1unif_vec = rep(NA, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)
logJointEvid0unif_vec = rep(NA, numSims)
logJointEvid1unif_vec = rep(NA, numSims)

# run simulations
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  # get x_uniform
  x_uniform_ind = sample(1:numx, N2)
  x_uniform_ind_mat[ , i] = x_uniform_ind
  x_uniform = x_seq[x_uniform_ind]
  x_uniform_mat[ , i] = x_uniform
  
  # generate mmed
  # step 1:
  step1_mmed = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step1_newpts = step1_mmed$addD
  
  x_train_new = c(x_train, step1_mmed$addD)
  x_train_ind_new = c(x_train_ind, step1_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 2:
  step2_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step2_newpts = step2_mmed$addD
  
  x_train_new = c(x_train_new, step2_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step2_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 3:
  step3_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step3_newpts = step3_mmed$addD
  
  x_train_new = c(x_train_new, step3_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step3_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  
  mmed_gp_list[[i]] = list(step1_mmed, step2_mmed, step3_mmed)
  
  # mmed inputs' predictions and evaluations
  newpts = c(step1_newpts, step2_newpts, step3_newpts)
  newpts_mat[ , i] = newpts
  newpts_ind = c(step1_mmed$indices, step2_mmed$indices, step3_mmed$indices)
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for mmed
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for mmed
  logJointEvid0mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # space-filling inputs' predictions and evaluations
  # (note: y_train is different for each function, which is why we need this in the loop)
  newpts = x_spacefill
  truey = y_seq[x_spacefill_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for spacefilling
  logPredEvid0sf_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1sf_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for spacefilling
  logJointEvid0sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # uniform inputs' predictions and evaluations
  newpts = x_uniform
  truey = y_seq[x_uniform_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for uniform
  RSS0unif_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1unif_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for uniform
  logPredEvid0unif_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1unif_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for uniform
  logJointEvid0unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
}

train3sims = list("grid" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train,
                  "x_train_ind" = x_train_ind,
                  # save the 3 designs (post input pts) that we're comparing
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "x_uniform" = x_uniform_mat,
                  "x_uniform_ind" = x_uniform_ind_mat,
                  "mmed_gp_list" = mmed_gp_list,
                  "x_mmed" = newpts_mat,
                  "x_mmed_ind" = newpts_ind_mat,
                  "y_mmed" = truey_mat,
                  # save evaluations for each of them also
                  "RSS0mmed_vec" = RSS0mmed_vec,
                  "RSS1mmed_vec" = RSS1mmed_vec,
                  "RSS0sf_vec" = RSS0sf_vec,
                  "RSS1sf_vec" = RSS1sf_vec,
                  "RSS0unif_vec" = RSS0unif_vec,
                  "RSS1unif_vec" = RSS1unif_vec,
                  "RSS01mmed" = RSS0mmed_vec / RSS1mmed_vec,
                  "RSS01sf" = RSS0sf_vec / RSS1sf_vec,
                  "RSS01unif" = RSS0unif_vec / RSS1unif_vec,
                  #
                  "logPredEvid0mmed_vec" = logPredEvid0mmed_vec,
                  "logPredEvid1mmed_vec" = logPredEvid1mmed_vec,
                  "logPredEvid0sf_vec" = logPredEvid0sf_vec,
                  "logPredEvid1sf_vec" = logPredEvid1sf_vec,
                  "logPredEvid0unif_vec" = logPredEvid0unif_vec,
                  "logPredEvid1unif_vec" = logPredEvid1unif_vec,
                  "logPredEvid01mmed" = logPredEvid0mmed_vec - logPredEvid1mmed_vec,
                  "logPredEvid01sf" = logPredEvid0sf_vec - logPredEvid1sf_vec,
                  "logPredEvid01unif" = logPredEvid0unif_vec - logPredEvid1unif_vec,
                  #
                  "logJointEvid0mmed_vec" = logJointEvid0mmed_vec,
                  "logJointEvid1mmed_vec" = logJointEvid1mmed_vec,
                  "logJointEvid0sf_vec" = logJointEvid0sf_vec,
                  "logJointEvid1sf_vec" = logJointEvid1sf_vec,
                  "logJointEvid0unif_vec" = logJointEvid0unif_vec,
                  "logJointEvid1unif_vec" = logJointEvid1unif_vec,
                  "logJointEvid01mmed" = logJointEvid0mmed_vec - logJointEvid1mmed_vec,
                  "logJointEvid01sf" = logJointEvid0sf_vec - logJointEvid1sf_vec,
                  "logJointEvid01unif" = logJointEvid0unif_vec - logJointEvid1unif_vec)
saveRDS(train3sims, paste(output_home,"mvp_seq_train3sims.rds", sep = ""))


## INPUT SET 4 ##

# input set 4
x_train_mat = matrix(NA, N, numSims)
x_train_ind_mat = matrix(NA, N, numSims)

# space-filling for input set 4
x_spacefill_ind = floor(c(1, 1 + ((numx - 1)/(N2 - 1)) * 1:((numx - 1) / ((numx - 1)/(N2 - 1)))))
x_spacefill = x_seq[x_spacefill_ind]

# these get rewritten for each train set
# uniform
x_uniform_mat = matrix(NA, N2, numSims)
x_uniform_ind_mat = matrix(NA, N2, numSims)
# mmed output list (just in case)
mmed_gp_list = list()
# mmed points (x and y)
newpts_mat = matrix(NA, N2, numSims)
newpts_ind_mat = matrix(NA, N2, numSims)
truey_mat = matrix(NA, N2, numSims)
# RSS
RSS0mmed_vec = rep(NA, numSims)
RSS1mmed_vec = rep(NA, numSims)
RSS0sf_vec = rep(NA, numSims)
RSS1sf_vec = rep(NA, numSims)
RSS0unif_vec = rep(NA, numSims)
RSS1unif_vec = rep(NA, numSims)
# log predictive density
logPredEvid0mmed_vec = rep(NA, numSims)
logPredEvid1mmed_vec = rep(NA, numSims)
logPredEvid0sf_vec = rep(NA, numSims)
logPredEvid1sf_vec = rep(NA, numSims)
logPredEvid0unif_vec = rep(NA, numSims)
logPredEvid1unif_vec = rep(NA, numSims)
# log joint density
logJointEvid0mmed_vec = rep(NA, numSims)
logJointEvid1mmed_vec = rep(NA, numSims)
logJointEvid0sf_vec = rep(NA, numSims)
logJointEvid1sf_vec = rep(NA, numSims)
logJointEvid0unif_vec = rep(NA, numSims)
logJointEvid1unif_vec = rep(NA, numSims)

# run simulations
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  # get x_train
  x_train_and_uniform_inds = sample(1:numx, Ntotal)
  x_train_ind = x_train_and_uniform_inds[1:N]
  x_train_ind_mat[1 , i] = x_train_ind
  x_train = x_seq[x_train_ind]
  x_train_mat[1 , i] = x_train
  # get y_train
  y_seq = y_seq_mat[ , i]
  y_train = y_seq[x_train_ind]
  # get x_uniform
  x_uniform_ind = x_train_and_uniform_inds[(N+1):Ntotal]
  x_uniform_ind_mat[ , i] = x_uniform_ind
  x_uniform = x_seq[x_uniform_ind]
  x_uniform_mat[ , i] = x_uniform
  
  # generate mmed
  # step 1:
  step1_mmed = add_MED_ms_oneatatime_data_gp(x_train, y_train, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step1_newpts = step1_mmed$addD
  
  x_train_new = c(x_train, step1_mmed$addD)
  x_train_ind_new = c(x_train_ind, step1_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 2:
  step2_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step2_newpts = step2_mmed$addD
  
  x_train_new = c(x_train_new, step2_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step2_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  # step 3:
  step3_mmed = add_MED_ms_oneatatime_data_gp(x_train_new, y_train_new, type01, l01, var_e = 1, N2 = stepN, 
                                             k = k, p = p, xmin = xmin, xmax = xmax, 
                                             nugget = nugget, alpha = alpha, buffer = 0, candidates = x_seq)
  step3_newpts = step3_mmed$addD
  
  x_train_new = c(x_train_new, step3_mmed$addD)
  x_train_ind_new = c(x_train_ind_new, step3_mmed$indices)
  y_train_new = y_seq[x_train_ind_new]
  
  mmed_gp_list[[i]] = list(step1_mmed, step2_mmed, step3_mmed)
  
  # mmed inputs' predictions and evaluations
  newpts = c(step1_newpts, step2_newpts, step3_newpts)
  newpts_mat[ , i] = newpts
  newpts_ind = c(step1_mmed$indices, step2_mmed$indices, step3_mmed$indices)
  newpts_ind_mat[ , i] = newpts_ind
  truey = y_seq[newpts_ind]
  truey_mat[ , i] = truey
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for mmed
  RSS0mmed_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1mmed_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for mmed
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid0mmed_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for mmed
  logJointEvid0mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1mmed_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # space-filling inputs' predictions and evaluations
  # (note: y_train is different for each function, which is why we need this in the loop)
  newpts = x_spacefill
  truey = y_seq[x_spacefill_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for spacefilling
  RSS0sf_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1sf_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for spacefilling
  logPredEvid0sf_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1sf_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for spacefilling
  logJointEvid0sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1sf_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
  
  # uniform inputs' predictions and evaluations
  newpts = x_uniform
  truey = y_seq[x_uniform_ind]
  H0_pred = getPredDistrSeq(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
  H1_pred = getPredDistrSeq(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
  postpredmu0 = H0_pred$pred_mean
  postpredmu1 = H1_pred$pred_mean
  # compute RSS0 and RSS1 for uniform
  RSS0unif_vec[i] = sum((postpredmu0 - truey)^2)  
  RSS1unif_vec[i] = sum((postpredmu1 - truey)^2)
  # compute log evidences (pred) for uniform
  logPredEvid0unif_vec[i] = dmvnorm(truey, mean = H0_pred$pred_mean, sigma = H0_pred$pred_var, log = TRUE)
  logPredEvid1unif_vec[i] = dmvnorm(truey, mean = H1_pred$pred_mean, sigma = H1_pred$pred_var, log = TRUE)
  # compute log evidences (joint) for uniform
  logJointEvid0unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[1], l01[1])
  logJointEvid1unif_vec[i] = logjointlik(newpts, truey, x_train, y_train, type01[2], l01[2])
}

RSS01mmed_vec = RSS0mmed_vec/RSS1mmed_vec
RSS01sf_vec = RSS0sf_vec/RSS1sf_vec

train4sims = list("grid" = x_seq,
                  "sim_fns" = y_seq_mat,
                  "x_train" = x_train_mat,
                  "x_train_ind" = x_train_ind_mat,
                  # save the 3 designs (post input pts) that we're comparing
                  "x_spacefill" = x_spacefill,
                  "x_spacefill_ind" = x_spacefill_ind,
                  "x_uniform" = x_uniform_mat,
                  "x_uniform_ind" = x_uniform_ind_mat,
                  "mmed_gp_list" = mmed_gp_list,
                  "x_mmed" = newpts_mat,
                  "x_mmed_ind" = newpts_ind_mat,
                  "y_mmed" = truey_mat,
                  # save evaluations for each of them also
                  "RSS0mmed_vec" = RSS0mmed_vec,
                  "RSS1mmed_vec" = RSS1mmed_vec,
                  "RSS0sf_vec" = RSS0sf_vec,
                  "RSS1sf_vec" = RSS1sf_vec,
                  "RSS0unif_vec" = RSS0unif_vec,
                  "RSS1unif_vec" = RSS1unif_vec,
                  "RSS01mmed" = RSS0mmed_vec / RSS1mmed_vec,
                  "RSS01sf" = RSS0sf_vec / RSS1sf_vec,
                  "RSS01unif" = RSS0unif_vec / RSS1unif_vec,
                  #
                  "logPredEvid0mmed_vec" = logPredEvid0mmed_vec,
                  "logPredEvid1mmed_vec" = logPredEvid1mmed_vec,
                  "logPredEvid0sf_vec" = logPredEvid0sf_vec,
                  "logPredEvid1sf_vec" = logPredEvid1sf_vec,
                  "logPredEvid0unif_vec" = logPredEvid0unif_vec,
                  "logPredEvid1unif_vec" = logPredEvid1unif_vec,
                  "logPredEvid01mmed" = logPredEvid0mmed_vec - logPredEvid1mmed_vec,
                  "logPredEvid01sf" = logPredEvid0sf_vec - logPredEvid1sf_vec,
                  "logPredEvid01unif" = logPredEvid0unif_vec - logPredEvid1unif_vec,
                  #
                  "logJointEvid0mmed_vec" = logJointEvid0mmed_vec,
                  "logJointEvid1mmed_vec" = logJointEvid1mmed_vec,
                  "logJointEvid0sf_vec" = logJointEvid0sf_vec,
                  "logJointEvid1sf_vec" = logJointEvid1sf_vec,
                  "logJointEvid0unif_vec" = logJointEvid0unif_vec,
                  "logJointEvid1unif_vec" = logJointEvid1unif_vec,
                  "logJointEvid01mmed" = logJointEvid0mmed_vec - logJointEvid1mmed_vec,
                  "logJointEvid01sf" = logJointEvid0sf_vec - logJointEvid1sf_vec,
                  "logJointEvid01unif" = logJointEvid0unif_vec - logJointEvid1unif_vec)
saveRDS(train4sims, paste(output_home,"mvp_seq_train4sims.rds", sep = ""))





