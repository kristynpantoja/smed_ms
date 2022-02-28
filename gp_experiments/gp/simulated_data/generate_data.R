for(scenario in c(1, 2, 4, 5)){
  ################################################################################
  # last updated: 2/27/2022
  # purpose: to test seqmedgp for scenarios 1 or 2
  #   where H1 is true
  
  ################################################################################
  # Sources/Libraries
  ################################################################################
  output_home = paste0("gp_experiments/gp/simulated_data/outputs")
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
  
  ################################################################################
  # simulation settings, shared for both scenarios
  ################################################################################
  
  # simulations settings
  numSims = 100
  xmin = -1
  xmax = 1
  numx = 10^3 + 1
  x_seq = seq(from = xmin, to = xmax, length.out = numx)
  sigmasq_measuremt = 1e-10
  sigmasq_signal = 1
  
  ################################################################################
  # Scenario settings
  ################################################################################
  
  if(scenario %in% c(1, 3, 6)){
    typeT = "matern"
    lT = 0.01
    pT = NULL
  } else if(scenario == 2){
    typeT = "periodic"
    lT = 0.01
    pT = 0.26
  } else if(scenario == 4){
    typeT = "periodic"
    lT = 0.5
    pT = 0.05
  } else if(scenario == 5){
    typeT = "squaredexponential"
    lT = 0.01
    pT = NULL
  } else{
    stop("invalid scenario given")
  }
  
  ################################################################################
  # generate functions 
  registerDoRNG(rng.seed)
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
  
  # plot GPs with space-filling design
  
  # # initial design
  # x_input_idx = ceiling(numx / 2)
  # x_input = x_seq[x_input_idx]
  # 
  # # choose the simulation and the design
  # idx = 1
  # y_seq = y_seq_mat[, idx]
  # y_input = y_seq[x_input_idx]
  # x_new_idx = sample(1:numx, 9, replace = FALSE)
  # x_new = x_seq[x_new_idx]
  # y_new = y_seq[x_new_idx]
  # 
  # # fit the true GP kernel
  # HT_predfn = getGPPredictive(x_seq, x_input, y_input, typeT, lT, 1, 
  #                             measurement.var = sigmasq_measuremt)
  # err = 2 * sqrt(diag(HT_predfn$pred_var))
  # ggdata = data.table(
  #   x = x_seq, 
  #   `True Function` = y_seq, 
  #   `HT Pred Mean` = HT_predfn$pred_mean, 
  #   lower = HT_predfn$pred_mean - err, 
  #   upper = HT_predfn$pred_mean + err
  # )
  # ggdata_pts = data.table(
  #   x = c(x_input, x_new), y = c(y_input, y_new), 
  #   color = c(rep(gg_color_hue(5)[2], length(x_input)), 
  #             rep(gg_color_hue(5)[1], length(x_new))), 
  #   shape = c(rep(8, length(x_input)), 
  #             rep(16, length(x_new)))
  # )
  # 
  # # plot it
  # plt = ggplot(ggdata) + 
  #   geom_path(aes(x = x, y = `True Function`)) + 
  #   geom_ribbon(aes(x = x, ymin = lower, ymax = upper), 
  #               color = gg_color_hue(5)[3], linetype = 2, 
  #               fill = gg_color_hue(5)[3], alpha = 0.1) +
  #   geom_path(aes(x = x, y = `HT Pred Mean`), color = gg_color_hue(5)[3]) +
  #   geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
  #              inherit.aes = FALSE, color = ggdata_pts$color, 
  #              shape = ggdata_pts$shape, 
  #              size = 3) +
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
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  } else{
    simulated_data_file = paste0(
      output_home,
      "/", typeT,
      "_l", lT,
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
