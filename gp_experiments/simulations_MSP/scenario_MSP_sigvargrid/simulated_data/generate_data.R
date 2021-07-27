################################################################################
# last updated: 05/27/2021
# purpose: to test seqmedgp for scenarios 1 or 2
#   where H1 is true

typeT = "periodic"
pT = 0.05
lT = 0.5
sigmasq_signal_seq = seq(0.5, 1.5, length.out = 101)
# https://towardsdatascience.com/understanding-gaussian-process-the-socratic-way-ba02369d804

################################################################################
# Sources/Libraries
################################################################################
output_home = paste0("gp_experiments/simulations_MSP/scenario_MSP_sigvargrid/simulated_data")
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMEDgp.R", sep = ""))
source(paste(functions_home, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_home, "/boxhill.R", sep = ""))
source(paste(functions_home, "/boxhill_gp.R", sep = ""))
source(paste(functions_home, "/kl_divergence.R", sep = ""))

library(mvtnorm)

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

library(ggplot2)
library(reshape2)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
library(data.table)

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 25
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10

################################################################################
# generate functions 
registerDoRNG(rng.seed)
foreach(
  i = 1:length(sigmasq_signal_seq)
) %dorng% {
  pT = 0.05
  lT = 0.5
  sigmasq_signal = sigmasq_signal_seq[i]
  null_cov = getCov(x_seq, x_seq, typeT, lT, pT, sigmasq_signal)
  null_mean = rep(0, numx)
  
  # the function values
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) 
  filename_append = ""
  if(is.null(sigmasq_measuremt)){
    y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov))
    filename_append = ""
  } else{
    y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov + 
                            sigmasq_measuremt * diag(numx))) 
    filename_append = "_noise"
  }
  
  # # choose the simulation to plot
  # idx = 1
  # # design pts (random, for illustration)
  # x_input_idx = sample(1:numx, 2, replace = FALSE)
  # x_input = x_seq[x_input_idx]
  # # data
  # y_seq = y_seq_mat[, idx]
  # y_input = y_seq[x_input_idx]
  # 
  # # fit the true GP kernel
  # HT_predfn = getGPPredictive(
  #   x_seq, x_input, y_input, typeT, lT, pT, 2, 
  #   measurement.var = sigmasq_measuremt)
  # err = 2 * sqrt(diag(HT_predfn$pred_var))
  # ggdata = data.table(
  #   x = x_seq, 
  #   `True Function` = y_seq, 
  #   `HT Pred Mean` = HT_predfn$pred_mean, 
  #   lower = HT_predfn$pred_mean - err, 
  #   upper = HT_predfn$pred_mean + err
  # )
  # 
  # # plot it
  # plt = ggplot(ggdata) + 
  #   geom_path(aes(x = x, y = `True Function`)) + 
  #   geom_ribbon(aes(x = x, ymin = lower, ymax = upper), 
  #               color = gg_color_hue(5)[3], linetype = 2, 
  #               fill = gg_color_hue(5)[3], alpha = 0.1) +
  #   geom_path(aes(x = x, y = `HT Pred Mean`), color = gg_color_hue(5)[3]) +
  #   geom_point(data = data.frame(
  #     x = x_input, y = y_input
  #   ), mapping = aes(x = x, y = y), 
  #   inherit.aes = FALSE, color = 2,
  #   size = 3) +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   labs(y = "y", x = "x", fill = "Function", color = "Function")
  # if(typeT == "periodic"){
  #   if(is.null(pT)) pT = 0.26
  #   plt = plt + geom_vline(xintercept = pT * (0:floor((xmax - xmin) / pT)), 
  #                          color = "blue", alpha = 0.5)
  # }
  # plot(plt)
  
  # save the results!
  if(typeT == "periodic"){
    if(is.null(pT)) pT = 0.26
    simulated_data_file = paste0(
      output_home,
      "/", typeT,
      "_l", lT,
      "_p", pT,
      "_sig", sigmasq_signal,
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  } else{
    simulated_data_file = paste0(
      output_home,
      "/", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  }
  saveRDS(
    list(
      x = x_seq, 
      null_mean = null_mean, 
      null_cov = null_cov, 
      numSims = numSims, 
      function_values_mat = y_seq_mat
    ), 
    file = simulated_data_file)
}
